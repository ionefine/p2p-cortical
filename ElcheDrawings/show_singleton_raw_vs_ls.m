function show_singleton_raw_vs_ls(vbl, d, r)
% Show the target, raw sim_draw subimage, and LS-fitted version
% for a singleton combo row r.
%
% Usage:
%   show_singleton_raw_vs_ls(vbl, 1, 1171)

if nargin < 2, d = 1; end
warning('off','MATLAB:rankDeficientMatrix');

% Paths and loads
drawDir   = fullfile(vbl.datadir, vbl.dirList(d).name);
modelsDir = fullfile(drawDir, 'models');
combosDir = fullfile(drawDir, 'combos');

if vbl.debugflag
    simFile   = fullfile(modelsDir,  [vbl.dirList(d).name, vbl.fileidstr, '_debug.mat']);
    comboFile = fullfile(combosDir,  [vbl.dirList(d).name, vbl.fileidstr, 'combo_debug.mat']);
else
    simFile   = fullfile(modelsDir,  [vbl.dirList(d).name, vbl.fileidstr, '.mat']);
    comboFile = fullfile(combosDir,  [vbl.dirList(d).name, vbl.fileidstr, 'combo.mat']);
end
S = load(simFile, 'sim_draw');   sim_draw = S.sim_draw;
C = load(comboFile, 'combo');    combo    = C.combo;

% Check row r is singleton
mask = combo.cmbx(r,:) > 0;
assert(sum(mask)==1, 'Row r must be a singleton combo row.');

% Build Xfull in same order as combine_sim_draws
ct = 1;
nPix  = size(sim_draw.subimg{1},1);
nCols = sim_draw.n_SubImages * sim_draw.nToSave;
Xfull = zeros(nPix, nCols, 'double');
for sd = 1:sim_draw.n_SubImages
    for ii = 1:sim_draw.nToSave
        Xfull(:, ct) = double(sim_draw.subimg{sd}(:, ii));
        ct = ct + 1;
    end
end

% Target and scaling
y = double(sim_draw.patient_img{1}(:));
scale = max(y) > 1.5;    % if 0..255, scale everything to 0..1
if scale
    y     = y / 255;
    Xfull = Xfull / 255;
end

% Resolve indices
j = find(mask, 1, 'first');                         % global column
sd_idx = ceil(j / sim_draw.nToSave);
ii_idx = j - (sd_idx-1) * sim_draw.nToSave;

% Raw subimage (as saved in sim_draw), before LS fit
raw_vec = double(sim_draw.subimg{sd_idx}(:, ii_idx));
if scale, raw_vec = raw_vec / 255; end

% LS-fitted version (what combo uses for singleton)
Xj = Xfull(:, j);
Bj = Xj \ y;
y_ls = Xj * Bj;

% Correlations
r_saved = sim_draw.corr{sd_idx}(ii_idx);   % raw corr stored earlier
r_raw   = corr_safe(raw_vec, y);           % recomputed raw corr (no fit)
r_ls    = corr_safe(y_ls,   y);            % LS corr (matches combo for singleton)

% Image shape
if ~isfield(sim_draw, 'size') || prod(sim_draw.size) ~= numel(y)
    L = numel(y); n = round(sqrt(L));
    assert(n*n == L, 'sim_draw.size missing or inconsistent.');
    sim_draw.size = [n n];
end
to2D = @(v) reshape(v, sim_draw.size);
norm01 = @(x) (x - min(x(:))) ./ max(1e-9, max(x(:)) - min(x(:)));

% Display
figure('Name','Target vs raw vs LS vs diff','Color','w');
tiledlayout(1,4,'Padding','compact','TileSpacing','compact');
nexttile; imshow(norm01(to2D(y)),    []); title('Target');
nexttile; imshow(norm01(to2D(raw_vec)), []); title(sprintf('Raw (sim\\_draw)\nr=%.3f (saved=%.3f)', r_raw, r_saved));
nexttile; imshow(norm01(to2D(y_ls)), []); title(sprintf('LS fit (combo)\nr=%.3f', r_ls));
nexttile; imshow(norm01(to2D(y_ls - raw_vec)), []); title('LS minus Raw');

fprintf('\nRow %d (global j=%d -> sd=%d, ii=%d)\n', r, j, sd_idx, ii_idx);
fprintf('corr_saved (sim_draw) : %.6f\n', r_saved);
fprintf('corr_raw   (recomp)   : %.6f\n', r_raw);
fprintf('corr_LS    (recomp)   : %.6f\n', r_ls);

end

function r = corr_safe(a,b)
a = a(:); b = b(:);
if all(a==a(1)) || all(b==b(1)), r = NaN; else, r = corr(a,b); end
end
