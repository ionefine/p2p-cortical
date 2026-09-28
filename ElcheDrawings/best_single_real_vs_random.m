function best_single_real_vs_random(vbl, d)

if nargin < 2, d = 1; end
warning('off','MATLAB:rankDeficientMatrix');

drawDir = fullfile(vbl.datadir, vbl.dirList(d).name);
modelsDir = fullfile(drawDir, 'models');
randDir   = fullfile(drawDir, 'random_models');

% load real
if vbl.debugflag
    simFile = fullfile(modelsDir, [vbl.dirList(d).name, vbl.fileidstr, '_debug.mat']);
else
    simFile = fullfile(modelsDir, [vbl.dirList(d).name, vbl.fileidstr, '2.mat']);
end
S = load(simFile, 'sim_draw');
sim_draw = S.sim_draw;

target_vec = double(sim_draw.patient_img{1}(:));

% determine shape
imgSizeKnown = isfield(sim_draw,'size') && numel(sim_draw.size)==2 ...
    && prod(sim_draw.size)==numel(target_vec);
if ~imgSizeKnown
    L = numel(target_vec);
    n = round(sqrt(L));
    if n*n ~= L
        error('Cannot infer image shape. Provide sim_draw.size in .mat.');
    end
    sim_draw.size = [n n];
end

normalize01 = @(x) (x - min(x(:))) ./ max(1e-9, (max(x(:)) - min(x(:))));

% target image
targetImg2D = reshape(target_vec, sim_draw.size);
if max(targetImg2D(:)) > 1.5, targetImg2D = targetImg2D/255; end
targetImg2D = normalize01(targetImg2D);

% best real
bestRealCorr = -inf; bestRealImg = []; bestRealInfo = struct('sd',[],'ii',[]);
for sd = 1:sim_draw.n_SubImages
    cvec = sim_draw.corr{sd};
    [cmax, ii] = max(cvec);
    if cmax > bestRealCorr
        bestRealCorr = cmax;
        bestRealImg  = double(sim_draw.subimg{sd}(:, ii));
        bestRealInfo.sd = sd; bestRealInfo.ii = ii;
    end
end
bestRealImg2D = reshape(bestRealImg, sim_draw.size);
if max(bestRealImg2D(:)) > 1.5, bestRealImg2D = bestRealImg2D/255; end
bestRealImg2D = normalize01(bestRealImg2D);

% -best rand per rep
bestRandPerRepCorr = nan(1, vbl.n_Reps);
bestRandPerRepImg  = cell(1, vbl.n_Reps);

for rep = 1:vbl.n_Reps
    try
        if vbl.debugflag
            rfile = fullfile(randDir, [vbl.dirList(d).name, vbl.fileidstr, '_', num2str(rep), '_rand_debug.mat']);
        else
            rfile = fullfile(randDir, [vbl.dirList(d).name, vbl.fileidstr, '_', num2str(rep), '_rand2.mat']);
        end
        if ~exist(rfile, 'file')
            fprintf('Skipping rep %d (file not found)\n', rep);
            continue
        end
        R = load(rfile, 'rand_draw');
        rand_draw = R.rand_draw;

        bestCorr = -inf; bestImg = [];
        for sd = 1:numel(rand_draw.subimg)
            cvec = rand_draw.corr{sd};
            [cmax, ii] = max(cvec);
            if cmax > bestCorr
                bestCorr = cmax;
                bestImg  = double(rand_draw.subimg{sd}(:, ii));
            end
        end

        if ~isempty(bestImg)
            img2D = reshape(bestImg, sim_draw.size);
            if max(img2D(:)) > 1.5, img2D = img2D/255; end
            bestRandPerRepImg{rep}  = normalize01(img2D);
            bestRandPerRepCorr(rep) = bestCorr;
        end
    catch ME
        fprintf('Rep %d load error: %s\n', rep, ME.message);
        continue
    end
end

% best random overall
validIdx = find(~isnan(bestRandPerRepCorr));
if isempty(validIdx)
    error('No valid random reps found.');
end
[bestRandCorr, bestIdx] = max(bestRandPerRepCorr(validIdx));
bestRandCorrRep = validIdx(bestIdx);
bestRandImg2D   = bestRandPerRepImg{bestRandCorrRep};

% display best images side by side
figure('Name','Target vs Top real vs Top random','Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% target
nexttile;
imshow(targetImg2D, []);
title('TARGET');

% real
nexttile;
imshow(bestRealImg2D, []);
title(sprintf('Top REAL (corr=%.3f)', bestRealCorr));

% random
nexttile;
imshow(bestRandImg2D, []);
title(sprintf('Top RANDOM rep=%d (corr=%.3f)', bestRandCorrRep, bestRandCorr));

% plot histogram
figure('Name','Best random per rep vs real','Color','w');
validCorrs = bestRandPerRepCorr(validIdx);
histogram(validCorrs, 'BinMethod','fd'); hold on;
xl = xline(bestRealCorr, 'r', 'LineWidth', 2);
xl.Label = sprintf('Real best = %.3f', bestRealCorr);
xl.LabelHorizontalAlignment = 'left';
xlabel('Correlation of best RANDOM (per rep)'); ylabel('Count');
title(sprintf('Random-best per rep vs Real-best (N reps = %d)', numel(validCorrs)));
grid on;

% final stats
fprintf('\nSummary for drawing %d:\n', d);
fprintf('  Best REAL corr: %.4f  (sd=%d, ii=%d)\n', bestRealCorr, bestRealInfo.sd, bestRealInfo.ii);
fprintf('  Best RANDOM corr overall: %.4f (rep=%d)\n', bestRandCorr, bestRandCorrRep);
fprintf('  Median best RANDOM corr across reps: %.4f\n\n', median(validCorrs));
end

