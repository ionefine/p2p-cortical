

function compare_real_vs_rand(vbl, d)
    % Load files
    draw_name = vbl.dirList(d).name;
    if vbl.debugflag
        combo_file = fullfile(vbl.datadir, draw_name, 'combos', ...
            [draw_name, vbl.fileidstr, '.mat']);
        rand_file  = fullfile(vbl.datadir, draw_name, 'combos', ...
            [draw_name, vbl.fileidstr, '.mat']);
    else
        combo_file = fullfile(vbl.datadir, draw_name, 'combos', ...
            [draw_name, vbl.fileidstr, '.mat']);
        rand_file  = fullfile(vbl.datadir, draw_name, 'combos', ...
            [draw_name, vbl.fileidstr, '_rand.mat']);
    end

    load(combo_file, 'combo');
    load(rand_file, 'rand_combo');

    % Look across nimg = 1,2,3
    for nimg = 1:3
        % Best real correlation for this nimg
        mask_real = (combo.nimg == nimg);
        best_real_corr = NaN;
        if any(mask_real)
            best_real_corr = max(combo.corr_val(mask_real));
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

        % Percentiles
        top5  = prctile(best_rand_per_rep, 95);
        top1  = prctile(best_rand_per_rep, 99);

        % Print summary
        fprintf('Image %s | nimg=%d\n', vbl.dirList(d).name, nimg);
        fprintf('Best REAL corr: %.4f\n', best_real_corr);
        fprintf('Random top 5%% threshold: %.4f\n', top5);
        fprintf('Random top 1%% threshold: %.4f\n', top1);

        % Compare where real sits
        if best_real_corr >= top1
            fprintf('real is within TOP 1%% of randoms.\n');
        elseif best_real_corr >= top5
            fprintf('real is within TOP 5%% of randoms.\n');
        else
            fprintf('real is BELOW top 5%% of randoms.\n');
        end
        fprintf('\n');
        percentile = sum(best_rand_per_rep <= best_real_corr) / numel(best_rand_per_rep ) * 100;
        fprintf('Real correlation is at the %.2f percentile\n', percentile);
    end
end