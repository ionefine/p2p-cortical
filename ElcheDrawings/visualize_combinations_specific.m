function visualize_combinations_specific(vbl, d)
tag = ed.iff(vbl.debugflag, '_debug', '');
simFile = fullfile(vbl.datadir, vbl.dirList(d).name, 'models', sprintf('%s%s%s.mat', vbl.dirList(d).name, vbl.fileidstr, tag));
cmbFile = fullfile(vbl.datadir, vbl.dirList(d).name, 'combos', sprintf('%s%s%s.mat', vbl.dirList(d).name, vbl.fileidstr, tag));

if ~isfile(simFile) || ~isfile(cmbFile)
    warning('Missing models or combos for %s', vbl.dirList(d).name);
    return;
end

S = load(simFile, 'sim_draw'); sim_draw = S.sim_draw;
C = load(cmbFile, 'combo');    combo = C.combo;

figure(1); clf; set(gcf,'Name', vbl.dirList(d).name);
for i = 1:numel(sim_draw.subimg)
    subplot(2,2,i);
    imagesc(sim_draw.patient_img{i});
    colormap(gray); axis off;
end

% Reconstruct on the fly (since combo.subimg omitted)
tmp_model = [];
ct = 1;
for sd = 1:sim_draw.n_SubImages
    for ii = 1:sim_draw.nToSave
        tmp_model.subimg(:, ct) = double(sim_draw.subimg{sd}(:, ii));
        ct = ct + 1;
    end
end
targetImg = single(sim_draw.patient_img{1});
target    = targetImg(:);

figure(20); clf;

% ----- best 3 first -----
idx3 = find(combo.nimg == 3);
[~, order3] = sort(combo.corr_val(idx3), 'descend');
best3_idx = idx3(order3(1));

mask3 = logical(combo.cmbx(best3_idx, :));
best3_src_idx = find(mask3);

X3 = single(tmp_model.subimg(:, mask3)) / 255;
B3 = X3 \ target;
recon3 = reshape(X3 * B3, size(targetImg));

% ----- best 1 from the best 3 components -----
best_corr1 = -Inf;
for a = 1:3
    src = best3_src_idx(a);
    X = single(tmp_model.subimg(:, src)) / 255;
    B = X \ target;
    recon = reshape(X * B, size(targetImg));
    r = corr(recon(:), double(target(:)));

    if r > best_corr1
        best_corr1 = r;
        recon1 = recon;
        best1_src_idx = src;
    end
end

% ----- best 2 from the best 3 components -----
pairs = nchoosek(best3_src_idx, 2);

best_corr2 = -Inf;
for p = 1:size(pairs,1)
    src = pairs(p,:);
    X = single(tmp_model.subimg(:, src)) / 255;
    B = X \ target;
    recon = reshape(X * B, size(targetImg));
    r = corr(recon(:), double(target(:)));

    if r > best_corr2
        best_corr2 = r;
        recon2 = recon;
        best2_src_idx = src;
    end
end

% ----- plot reconstructions -----
subplot(2,3,1);
imagesc(recon1); axis image off; colormap gray;
title(sprintf('1 of best 3\nidx=%d, r=%.3f', best1_src_idx, best_corr1));

subplot(2,3,2);
imagesc(recon2); axis image off; colormap gray;
title(sprintf('2 of best 3\nidx=%s, r=%.3f', mat2str(best2_src_idx), best_corr2));

subplot(2,3,3);
imagesc(recon3); axis image off; colormap gray;
title(sprintf('Best 3\nidx=%s, r=%.3f', mat2str(best3_src_idx), combo.corr_val(best3_idx)));

% ----- plot individual components from best 3 -----
for j = 1:3
    img = reshape(tmp_model.subimg(:, best3_src_idx(j)), size(targetImg));

    subplot(2,3,3+j);
    imagesc(img); axis image off; colormap gray;
    title(sprintf('x_%d idx=%d', j, best3_src_idx(j)));
end
end