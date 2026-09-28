function sig_level = compare_real_vs_rand2(vbl, d)

% Load files
draw_name = vbl.dirList(d).name;

if vbl.debugflag
    combo_file = fullfile(vbl.datadir, draw_name, 'combos', ...
        [draw_name, vbl.fileidstr, '.mat']);

    rand_file = fullfile(vbl.datadir, draw_name, 'combos', ...
        [draw_name, vbl.fileidstr, '.mat']);
else
    combo_file = fullfile(vbl.datadir, draw_name, 'combos', ...
        [draw_name, vbl.fileidstr, '.mat']);

    rand_file = fullfile(vbl.datadir, draw_name, 'combos', ...
        [draw_name, vbl.fileidstr, '_rand.mat']);
end

% Skip if files missing
if ~isfile(combo_file)
    fprintf('Missing combo file for %s. Skipping.\n', draw_name);
    sig_level = NaN;
    return
end

if ~isfile(rand_file)
    fprintf('Missing random file for %s. Skipping.\n', draw_name);
    sig_level = NaN;
    return
end

% Load real combo
S1 = load(combo_file);
if ~isfield(S1, 'combo')
    fprintf('No combo variable in %s. Skipping.\n', combo_file);
    sig_level = NaN;
    return
end
combo = S1.combo;

% Load random combo
S2 = load(rand_file);

if isfield(S2, 'rand_combo')
    rand_combo = S2.rand_combo;
elseif isfield(S2, 'combo')
    rand_combo = S2.combo;
else
    fprintf('No rand_combo or combo variable in %s. Skipping.\n', rand_file);
    sig_level = NaN;
    return
end

% Store results for nimg = 1,2,3
real_corr = NaN(3,1);
rand_top5 = NaN(3,1);
passes_95 = false(3,1);

for nimg = 1:3

    % Best real correlation for this nimg
    mask_real = (combo.nimg == nimg);
    if any(mask_real)
        real_corr(nimg) = max(combo.corr_val(mask_real));
    end

    % Best random correlation per rep
    best_rand_per_rep = NaN(vbl.n_Reps,1);

    for rep = 1:vbl.n_Reps
        mask_rand = (rand_combo.nimg(rep,:) == nimg);

        if any(mask_rand)
            best_rand_per_rep(rep) = max(rand_combo.corr_val(rep, mask_rand));
        end
    end

    best_rand_per_rep = best_rand_per_rep(~isnan(best_rand_per_rep));

    if ~isempty(best_rand_per_rep)
        rand_top5(nimg) = prctile(best_rand_per_rep, 95);
    end

    % Did real beat the 95th percentile?
    if ~isnan(real_corr(nimg)) && ~isnan(rand_top5(nimg))
        passes_95(nimg) = real_corr(nimg) >= rand_top5(nimg);
    end
end

% Assign highest level: prefer 3, then 2, otherwise 1
if passes_95(3)
    sig_level = 3;
elseif passes_95(2)
    sig_level = 2;
else
    sig_level = 1;
end

fprintf('Image %s assigned sig_level = %d\n', draw_name, sig_level);
fprintf('Passes 95th percentile: nimg1=%d, nimg2=%d, nimg3=%d\n', ...
    passes_95(1), passes_95(2), passes_95(3));

% Save result
save_file = fullfile(vbl.datadir, draw_name, 'combos', ...
    [draw_name, vbl.fileidstr, '_real_vs_rand_summary.mat']);

save(save_file, 'draw_name', 'real_corr', 'rand_top5', 'passes_95', 'sig_level');

end