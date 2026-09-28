function audit_singleton_combos(vbl, d, tol, showExample)
% Verify that 1-image combos in combo.mat match original single-image stats.
% Compares three numbers for each singleton:
%   1) combo.corr_val on the singleton row (from combine_sim_draws)
%   2) recomputed corr using LS with that single column
%   3) sim_draw.corr{sd}(ii) stored earlier
%
% Usage:
%   audit_singleton_combos(vbl, 1);
%   audit_singleton_combos(vbl, 1, 1e-8, true);

if nargin < 2, d = 1; end
if nargin < 3, tol = 1e-8; end
if nargin < 4, showExample = true; end
warning('off','MATLAB:rankDeficientMatrix');

% Paths
drawDir   = fullfile(vbl.datadir, vbl.dirList(d).name);
modelsDir = fullfile(drawDir, 'models');
combosDir = fullfile(drawDir, 'combos');

% Load sim_draw and combo
if vbl.debugflag
    simFile   = fullfile(modelsDir,  [vbl.dirList(d).name, vbl.fileidstr, '_debug.mat']);
    comboFile = fullfile(combosDir,  [vbl.dirList(d).name, vbl.fileidstr, 'combo_debug.mat']);
else
    simFile   = fullfile(modelsDir,  [vbl.dirList(d).name, vbl.fileidstr, '2.mat']);
    comboFile = fullfile(combosDir,  [vbl.dirList(d).name, vbl.fileidstr, 'combo2.mat']);
end
S = load(simFile, 'sim_draw');   sim_draw = S.sim_draw;
C = load(comboFile, 'combo');    combo    = C.combo;

% Build Xfull in the exact column order used for combos
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
if max(y) > 1.5
    y     = y / 255;
    Xfull = Xfull / 255;
end

% Find singleton rows
isSingleton = (sum(combo.cmbx, 2) == 1);
rows = find(isSingleton);
if isempty(rows)
    fprintf('No singleton combos found for drawing %d.\n', d);
    return
end
fprintf('Found %d singleton combo rows for drawing %d.\n', numel(rows), d);

% Helpers
normalize01 = @(x) (x - min(x(:))) ./ max(1e-9, (max(x(:)) - min(x(:))));
if ~isfield(sim_draw, 'size') || prod(sim_draw.size) ~= numel(y)
    L = numel(y); n = round(sqrt(L));
    assert(n*n == L, 'Cannot infer sim_draw.size; please save it in sim_draw.');
    sim_draw.size = [n n];
end

% Storage for comparisons
diff_combo_vs_recompute = zeros(numel(rows), 1);
diff_combo_vs_savedcorr = zeros(numel(rows), 1);
report = strings(numel(rows), 1);

% Check each singleton
for k = 1:numel(rows)
    r = rows(k);
    mask = combo.cmbx(r, :) > 0;
    j = find(mask, 1, 'first');           % global column index

    % Recompute solo corr via LS
    Xj   = Xfull(:, j);
    Bj   = Xj \ y;
    yhat = Xj * Bj;
    rc_recompute = safe_corr(yhat, y);

    % Combo stored corr for this singleton row
    rc_combo = combo.corr_val(r);

    % Map global j -> (sd, ii) to fetch sim_draw.corr{sd}(ii)
    [sd_idx, ii_idx] = global_to_sdii(j, sim_draw.n_SubImages, sim_draw.nToSave);
    rc_saved = sim_draw.corr{sd_idx}(ii_idx);

    % Record diffs
    diff_combo_vs_recompute(k) = abs(rc_combo - rc_recompute);
    diff_combo_vs_savedcorr(k) = abs(rc_combo - rc_saved);

    report(k) = sprintf('row=%d  j=%d  (sd=%d,ii=%d)  corr_combo=%.6f  corr_recomp=%.6f  corr_saved=%.6f  |diffs|=[%.2e, %.2e]', ...
        r, j, sd_idx, ii_idx, rc_combo, rc_recompute, rc_saved, ...
        diff_combo_vs_recompute(k), diff_combo_vs_savedcorr(k));
end

% Print summary
fprintf('\nSummary (tol = %.1e):\n', tol);
n_bad1 = sum(diff_combo_vs_recompute > tol);
n_bad2 = sum(diff_combo_vs_savedcorr > tol);
fprintf('  combo vs recomputed:  max |diff| = %.3e  (violations: %d of %d)\n', ...
    max(diff_combo_vs_recompute), n_bad1, numel(rows));
fprintf('  combo vs sim_draw:    max |diff| = %.3e  (violations: %d of %d)\n', ...
    max(diff_combo_vs_savedcorr), n_bad2, numel(rows));

% Print offending lines if any
if n_bad1 > 0 || n_bad2 > 0
    fprintf('\nMismatches beyond tolerance:\n');
    for k = 1:numel(rows)
        if diff_combo_vs_recompute(k) > tol || diff_combo_vs_savedcorr(k) > tol
            fprintf('%s\n', report(k));
        end
    end
else
    fprintf('All singleton combos match within tolerance.\n');
end

% Optional quick visual for the first singleton
if showExample
    r = rows(1);
    mask = combo.cmbx(r, :) > 0;
    j = find(mask, 1, 'first');
    Xj = Xfull(:, j);
    Bj = Xj \ y;
    yhat = Xj * Bj;

    figure('Name','Singleton audit: target vs solo', 'Color','w');
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    nexttile; imshow(normalize01(reshape(y,    sim_draw.size)), []); title('Target');
    nexttile; imshow(normalize01(reshape(yhat, sim_draw.size)), []); title('Singleton recon');
end

end

% ---- helpers ----
function r = safe_corr(a, b)
    a = a(:); b = b(:);
    if all(a == a(1)) || all(b == b(1)), r = NaN; else, r = corr(a, b); end
end

function [sd, ii] = global_to_sdii(j, n_SubImages, nToSave)
    sd = ceil(j / nToSave);
    ii = j - (sd - 1) * nToSave;
    sd = max(1, min(sd, n_SubImages));
    ii = max(1, min(ii, nToSave));
end
