% === inputs you set ===
d = 2;  % which drawing index

% === load matching files for THIS drawing ===
if vbl.debugflag
    S1 = load(fullfile(vbl.datadir, vbl.dirList(d).name, 'combos', ...
        [vbl.dirList(d).name, vbl.fileidstr, 'combo_debug.mat']), 'combo');
    S2 = load(fullfile(vbl.datadir, vbl.dirList(d).name, 'combos', ...
        [vbl.dirList(d).name, vbl.fileidstr, 'combo_rand_debug.mat']), 'rand_combo');
    S3 = load(fullfile(vbl.datadir, vbl.dirList(d).name, 'models', ...
        [vbl.dirList(d).name, vbl.fileidstr, '_debug.mat']), 'sim_draw');
else
    S1 = load(fullfile(vbl.datadir, vbl.dirList(d).name, 'combos', ...
        [vbl.dirList(d).name, vbl.fileidstr, 'combo.mat']), 'combo');
    S2 = load(fullfile(vbl.datadir, vbl.dirList(d).name, 'combos', ...
        [vbl.dirList(d).name, vbl.fileidstr, 'combo_rand.mat']), 'rand_combo');
    S3 = load(fullfile(vbl.datadir, vbl.dirList(d).name, 'models', ...
        [vbl.dirList(d).name, vbl.fileidstr, '.mat']), 'sim_draw');
end
combo = S1.combo;
rand_combo = S2.rand_combo;
sim_draw = S3.sim_draw;

% true image shape for reshaping
sz = size(sim_draw.patient_img{1});
target = double(sim_draw.patient_img{1}(:));

% build "real" subimage pool in same order as combine_sim_draws
realimg = [];
for sd = 1:sim_draw.n_SubImages
    for ii = 1:sim_draw.nToSave
        realimg = [realimg double(sim_draw.subimg{sd}(:,ii))/255]; %#ok<AGROW>
    end
end

figure('Color','w','Name',sprintf('Drawing %d best over reps', d));

for nimg_val = 1:3
    %% best non-random for this nimg
    idx_non = find(combo.nimg == nimg_val);
    if ~isempty(idx_non)
        [best_corr_real, best_idx_rel] = max(combo.corr_val(idx_non));
        real_vec = combo.subimg(idx_non(best_idx_rel), :);
        % reshape with the true patient image size
        real_best_img = reshape(real_vec, sz);
    else
        best_corr_real = NaN;
        real_best_img = nan(sz);
    end

    %% best random across ALL reps for this nimg
    mask_rand = (rand_combo.nimg == nimg_val);
    if any(mask_rand(:))
        % locate the absolute best corr across reps for this nimg size
        corr_mat = rand_combo.corr_val;
        corr_mat(~mask_rand) = -Inf;
        [best_corr_rand, linear_idx] = max(corr_mat(:));
        [best_rep, best_combo_idx] = ind2sub(size(corr_mat), linear_idx);

        % load the raw random draw for that repetition
        if vbl.debugflag
            R = load(fullfile(vbl.datadir, vbl.dirList(d).name, 'random_models', ...
                [vbl.dirList(d).name, vbl.fileidstr, '_', num2str(best_rep), '_rand_debug.mat']), 'rand_draw');
        else
            R = load(fullfile(vbl.datadir, vbl.dirList(d).name, 'random_models', ...
                [vbl.dirList(d).name, vbl.fileidstr, '_', num2str(best_rep), '_rand.mat']), 'rand_draw');
        end
        rand_draw = R.rand_draw;

        % flatten random subimages (note: this may not reproduce any per-rep shuffles)
        randimg = [];
        for sd = 1:numel(rand_draw.subimg)
            randimg = [randimg double(rand_draw.subimg{sd})/255]; %#ok<AGROW>
        end

        % reconstruct using the saved mask
        mask = logical(rand_combo.cmbx(best_combo_idx, :));
        Xr = randimg(:, mask);
        if sum(mask) > 1
            Xr(:,1) = realimg(:,1);  % anchor rule used in combine_random_models
        end
        B = Xr \ target;
        recon_rand = Xr * B;
        rand_best_img = reshape(recon_rand, sz);
    else
        best_corr_rand = NaN;
        rand_best_img = nan(sz);
    end

    %% plot side by side
    subplot(3,2,2*nimg_val-1);
    imagesc(real_best_img); axis image off; colormap gray
    title(sprintf('Real best (%d img) corr=%.3f', nimg_val, best_corr_real));

    subplot(3,2,2*nimg_val);
    imagesc(rand_best_img); axis image off; colormap gray
    title(sprintf('Rand best over reps (%d img) corr=%.3f', nimg_val, best_corr_rand));
end
