for d = 1; %:vbl.n_Drawings
    % load non-random and random combos
    if vbl.debugflag
        load(fullfile(vbl.datadir, vbl.dirList(d).name, 'combos', ...
            [vbl.dirList(d).name, vbl.fileidstr, 'combo_debug.mat']), 'combo');
        load(fullfile(vbl.datadir, vbl.dirList(d).name, 'combos', ...
            [vbl.dirList(d).name, vbl.fileidstr, 'combo_rand_debug.mat']), 'rand_combo');
    else
        load(fullfile(vbl.datadir, vbl.dirList(d).name, 'combos', ...
            [vbl.dirList(d).name, vbl.fileidstr, 'combo.mat']), 'combo');
        load(fullfile(vbl.datadir, vbl.dirList(d).name, 'combos', ...
            [vbl.dirList(d).name, vbl.fileidstr, 'combo_rand.mat']), 'rand_combo');
    end

    figure('Name', sprintf('Drawing %d', d), 'Color', 'w');
    nimg_vals = [1 2 3];  % only combos of 1, 2, 3

    for i = 1:length(nimg_vals)
        nimg_val = nimg_vals(i);

        % --- non-random best ---
        idx_nonrand = (combo.nimg == nimg_val);
        if any(idx_nonrand)
            best_corr_nonrand = max(combo.corr_val(idx_nonrand));
        else
            best_corr_nonrand = NaN;
        end

        % --- random best per rep ---
        best_corr_rand = NaN(vbl.n_Reps, 1);
        for rep = 1:vbl.n_Reps
            idx_rand = (rand_combo.nimg(rep,:) == nimg_val);
            if any(idx_rand)
                best_corr_rand(rep) = max(rand_combo.corr_val(rep, idx_rand));
            end
        end

        % --- plot histogram ---
subplot(1,3,i);

vals = best_corr_rand(~isnan(best_corr_rand));
all_vals = [vals(:); best_corr_nonrand];   % combine random + nonrandom

if ~isempty(all_vals)
    lo = min(all_vals);
    hi = max(all_vals);
    pad = 0.05 * (hi - lo + eps);          % 5% padding
    edges = linspace(lo - pad, hi + pad, 20);  % ~20 bins
    
    histogram(vals, 'BinEdges', edges, 'FaceColor',[0.2 0.6 0.8]); hold on;
    if ~isnan(best_corr_nonrand)
        xline(best_corr_nonrand, 'r-', 'LineWidth',2);
    end
    xlim([lo - pad, hi + pad]);
else
    histogram(0, 'BinEdges', [0 1]); xlim([0 1]);
end

title(sprintf('%d image(s)', nimg_val));
xlabel('Correlation');
ylabel('Count');

        
    end

    sgtitle(sprintf('Drawing %d: Top random vs non-random', d));
end