% ed.m
%
% Support functions for the Elche Drawings Project.
% Usage: call as static class methods, e.g., vbl = ed.setup('Z');
%
% Major groups:
% - Setup and I/O helpers
% - Patient drawing loading and cropping
% - Cortical/visual model definition and phosphene generation
% - Simulation of "best" singletons and saving
% - Random phosphene generation, combination, and analysis
% - Visualization utilities

classdef ed
    methods(Static)

        % ---------------------------
        % Setup
        % ---------------------------
        function vbl = setup(loc)
            % Setup paths and locate data directory
            % loc options:
            %   'Z'      -> hardcoded network drive example
            %   'auto'   -> use pwd
            %   PC paths -> pass a full path string
            %
            % Populates:
            %   vbl.datadir, vbl.dirList, vbl.n_Drawings, vbl.n_SimDraw,
            %   vbl.n_DrawCombined, vbl.debugflag, vbl.n_Reps, vbl.rectify,
            %   vbl.fileidstr, vbl.cleanstartflag

            % Add helper paths relative to this file
            here = fileparts(mfilename('fullpath'));
            addpath(here);
            addpath(fullfile(here, '..'));
            addpath(fullfile(here, '..', 'imregcode'));

            if nargin < 1
                loc = 'auto';
            end

            arch = computer('arch');
            switch loc
                case 'Z'
                    vbl.datadir = fullfile('Z:', 'IoneFine', 'MyProjects', 'p2p-cortical', 'ElcheDrawings');
                case 'auto'
                    vbl.datadir = pwd;
                otherwise
                    % If a full path is supplied, use it
                    if isfolder(loc)
                        vbl.datadir = loc;
                    else
                        % Fallback: on Linux/GLNXA64 or Windows, use pwd
                        vbl.datadir = pwd;
                    end
            end

            % Discover drawings
            vbl.dirList = dir(fullfile(vbl.datadir, '*img*'));
            vbl.n_Drawings = numel(vbl.dirList);
            if vbl.n_Drawings == 0
                error('No drawing folders found in %s (pattern "*img*").', vbl.datadir);
            end

            % Configuration
            vbl.n_SimDraw      = 24; % total singletons to keep across subimages
            vbl.n_DrawCombined = 3;  % max images combined in a set

            vbl.debugflag = 1;
            vbl.n_Reps    = iff(vbl.debugflag==1, 4, 100); % helper below

            vbl.rectify   = 1;   % 1: ignore negatives when scaling, 0: bipolar scaling
            vbl.fileidstr = '_10_2_2025';

            % Clean start (ask user unless overridden upstream)
            vbl.cleanstartflag = 1; % default ask
            vbl.cleanstartflag = input('Delete all previous simulations? 1=yes, 0=no ... ');
        end

        % ---------------------------
        % Patient drawings I/O
        % ---------------------------
        function p_draw = load_patient_drawings(vbl)
            % Load each drawing's cropped subimages and preallocate arrays

            for d = 1:vbl.n_Drawings
                drawDir = fullfile(vbl.datadir, vbl.dirList(d).name);
                disp(['loading ', drawDir]);

                % Optional: clean previous results inside this drawing folder
                if vbl.cleanstartflag
                    ed.ensure_dirs(drawDir); % creates needed subfolders
                    ed.clean_dir_safe(fullfile(drawDir, 'combos'));
                    ed.clean_dir_safe(fullfile(drawDir, 'models'));
                    ed.clean_dir_safe(fullfile(drawDir, 'random_models'));
                    ed.clean_dir_safe(fullfile(drawDir)); % *.mat at root
                else
                    ed.ensure_dirs(drawDir);
                end

                subImgs = dir(fullfile(drawDir, 'drawings', '*edit_crop_*.png'));
                p_draw(d).n_SubImages = numel(subImgs);
                if p_draw(d).n_SubImages == 0
                    error('No subimages found in %s/drawings', drawDir);
                end

                p_draw(d).nToSave = ceil(vbl.n_SimDraw / p_draw(d).n_SubImages);

                for sd = 1:p_draw(d).n_SubImages
                    fname = fullfile(drawDir, 'drawings', sprintf('%s_edit_crop_%d.png', vbl.dirList(d).name, sd-1));

                    tmp = importdata(fname);
  
                    if isstruct(tmp)
                        img = uint8(mean(tmp.cdata, 3));
                    else
                        img = uint8(mean(tmp, 3));
                    end

                p_draw(d).patient_img{sd} = img;
                p_draw(d).size{sd}        = size(img);
                if sd == 1
                    p_draw(d).ref2d = imref2d(size(img)); % reference for warping
                end

                % Preallocate top-K arrays (per subimage)
                p_draw(d).sim_img{sd} = uint8(NaN([size(img), vbl.n_SimDraw])); % for potential storage
                p_draw(d).corr{sd}    = -inf(p_draw(d).nToSave, 1);
                p_draw(d).ssim{sd}    = -inf(p_draw(d).nToSave, 1);
                p_draw(d).radius{sd}  = NaN(p_draw(d).nToSave, 1);
                p_draw(d).x{sd}       = NaN(p_draw(d).nToSave, 1);
                p_draw(d).y{sd}       = NaN(p_draw(d).nToSave, 1);
                p_draw(d).subID{sd}   = NaN(p_draw(d).nToSave, 1);
                p_draw(d).subimg{sd}  = NaN(numel(img), p_draw(d).nToSave); % store vectorized images
            end
        end
    end

    function crop_drawings(vbl)
        % Create cropped versions for each drawing if not already present.
        for dd = 1:vbl.n_Drawings
            drawName = vbl.dirList(dd).name;
            drawDir  = fullfile(vbl.datadir, drawName, 'drawings');

            origFile = fullfile(drawDir, [drawName, '.png']);
            editFile = fullfile(drawDir, [drawName, '_edit.png']);

            if ~isfile(origFile) || ~isfile(editFile)
                warning('Missing original or edited image in %s', drawDir);
                continue;
            end

            orig_img = double(imread(origFile));
            edit_img = double(imread(editFile));

            sz  = size(orig_img);
            img = mean(orig_img, 3);

            % Crop white space via zero rows/cols along midlines
            rowMid = round(sz(1)/2);
            colMid = round(sz(2)/2);

            cidx = find(img(rowMid, :) == 0);
            ridx = find(img(:, colMid) == 0);
            if numel(cidx) < 2 || numel(ridx) < 2
                warning('Failed to auto-crop %s; saving originals as crops.', drawName);
                orig_img_tmp = orig_img;
                edit_img_tmp = edit_img;
            else
                crop = [ridx(1)+1, cidx(1)+1, ridx(end)-1, cidx(end)-1]; % [r1 c1 r2 c2]
                orig_img_tmp = orig_img(crop(1):crop(3), crop(2):crop(4), :);
                edit_img_tmp = edit_img(crop(1):crop(3), crop(2):crop(4), :);
            end

            % If tighter crop function is not used, save the tmp
            orig_out = fullfile(drawDir, [drawName, '_crop.png']);
            edit_out = fullfile(drawDir, [drawName, '_edit_crop_0.png']);

            % Normalize to [0,1] for saving
            imwrite(ed.norm01(orig_img_tmp), orig_out);
            imwrite(ed.norm01(edit_img_tmp), edit_out);
        end
    end

    % ---------------------------
    % Cortical/visual model
    % ---------------------------
    function [c, v, trl, tp] = define_cortical_model(vbl)
        disp('Defining the cortical model...');

        tp = p2p_c.define_temporalparameters();

        trl.freq = NaN;              % decouple intensity from params
        trl      = p2p_c.define_trial(tp, trl);

        if vbl.debugflag
            v.pixperdeg = 4;
            c.pixpermm  = 4;
            nElectrodes = 20;
        else
            v.pixperdeg = 10;
            c.pixpermm  = 10;
            nElectrodes = 1000;
        end

        c.cortexHeight = [-30, 30];
        c.cortexLength = [-55, 55];
        c.onoff_ratio  = 0.9;

        v.visfieldHeight = [-30, 30];
        v.visfieldWidth  = [-30, 30];

        v = p2p_c.define_visualmap(v);
        c = p2p_c.define_cortex(c);
        [c, v] = p2p_c.generate_corticalmap(c, v);

        % Electrode placement
        theta = 2*pi*rand(1, nElectrodes);
        r     = exp(4*rand(1, nElectrodes)); r = r*3./max(r(:));
        ecc   = exp((log(15) * rand(1, nElectrodes)));

        c.I_k = 1000; % electric field falloff constant

        [x, y] = pol2cart(theta(:), ecc(:));
        iperm  = randperm(numel(x));

        for i = 1:numel(x)
            v.e(i).x = x(iperm(i));
            v.e(i).y = y(iperm(i));
            c.e(i).radius = r(iperm(i));
        end

        disp(['Defined ', num2str(numel(x)), ' electrodes.']);
        c = p2p_c.define_electrodes(c, v);
    end

    function p_draw = simulate_drawings(c_orig, v_orig, trl, tp, p_draw, vbl)
        % For each electrode: generate a phosphene, align to each subimage,
        % and keep top-K per subimage with diversity constraint.

        corrThresh = 0.90; % limit near-duplicate phosphenes

        for eIdx = 1:numel(c_orig.e)
            disp(sprintf('Electrode %d / %d', eIdx, numel(c_orig.e)));

            % Restrict to single electrode
            c = ed.safe_rmfield(c_orig, {'e'}); v = ed.safe_rmfield(v_orig, {'e'});
            c.e = c_orig.e(eIdx);
            v.e = v_orig.e(eIdx);

            c = p2p_c.generate_ef(c);
            v = p2p_c.generate_corticalelectricalresponse(c, v);

            img = uint8(ed.generate_phosphene(v, tp, trl, vbl)); % 2D uint8

            for d = 1:vbl.n_Drawings
                for sd = 1:p_draw(d).n_SubImages
                    % Align img to patient subimage
                    target = p_draw(d).patient_img{sd};
                    [s, r] = findScaleRotationNGC_if(single(img), single(target));
                    [tform, peakcorr] = resolveSimilarityRotationAmbiguityNGC_if(single(img), single(target), s, r);
                    img_aligned = imwarp(img, tform, "OutputView", p_draw(d).ref2d);

                    if peakcorr <= min(p_draw(d).corr{sd})
                        continue; % not promising vs current pool
                    end

                    % Normalize and vectorize
                    vImg = double(img_aligned(:));
                    vmax = max(vImg); vmin = min(vImg);
                    if vbl.rectify
                        vImg = max(vImg - vmin, 0);
                        if vmax > 0
                            vImg = vImg / vmax * 255;
                        end
                        scaled_vec = uint8(vImg);
                    else
                        denom = max(vmax - vmin, eps);
                        scaled_vec = int8(128 * (vImg - vmin) / denom);
                    end

                    % Dimensions check
                    if numel(scaled_vec) ~= numel(target)
                        error('Model and patient images are different sizes for d=%d sd=%d.', d, sd);
                    end

                    % Diversity check vs saved pool
                    if eIdx == 1 || all(isnan(p_draw(d).subimg{sd}(:)))
                        allcorr = 0;
                    else
                        pool = p_draw(d).subimg{sd}; % [Npix x K]
                        allcorr = zeros(1, size(pool,2));
                        for k = 1:size(pool,2)
                            if all(isnan(pool(:,k)))
                                allcorr(k) = 0;
                            else
                                allcorr(k) = corr(double(pool(:,k)), double(scaled_vec));
                            end
                        end
                    end

                    % Select replacement index
                    if max(allcorr) > corrThresh
                        [~, ridx] = max(allcorr);
                    else
                        [~, ridx] = min(p_draw(d).corr{sd});
                    end

                    % Compute metrics vs target
                    newCorr = corr(double(target(:)), double(scaled_vec(:)));
                    newSSIM = ssim(reshape(double(scaled_vec), size(target)), double(target), 'Exponents', [0 0 1]);

                    if newCorr > p_draw(d).corr{sd}(ridx)
                        p_draw(d).subimg{sd}(:, ridx) = double(scaled_vec(:));
                        p_draw(d).radius{sd}(ridx)   = c.e.radius;
                        p_draw(d).x{sd}(ridx)        = v.e.x;
                        p_draw(d).y{sd}(ridx)        = v.e.y;
                        p_draw(d).corr{sd}(ridx)     = newCorr;
                        p_draw(d).ssim{sd}(ridx)     = newSSIM;
                        p_draw(d).subID{sd}(ridx)    = sd;
                        p_draw(d).size{sd}           = size(img_aligned);
                    end
                end
            end
        end
    end

    function img = generate_phosphene(v, tp, trl, vbl)
        % Generate a single phosphene image for the current electrode in v
        trl = p2p_c.generate_phosphene(v, tp, trl);
        img = max(trl.max_phosphene, [], 3);
        img = crop_img(img, 20);
        img = img ./ max(abs(img(:)) + eps);

        if vbl.rectify
            img(img < 0) = 0;
            img = uint8(255 * img);
        else
            img = uint8(127 * (img + 0.5)); % bipolar scaling to ~[0,255]
        end
    end

    function save_simulated_drawings(p_draw, vbl)
        for d = 1:vbl.n_Drawings
            sim_draw.n_SubImages = p_draw(d).n_SubImages;
            sim_draw.nToSave     = p_draw(d).nToSave;
            sim_draw.patient_img = p_draw(d).patient_img;
            sim_draw.ref2d       = p_draw(d).ref2d;

            for sd = 1:sim_draw.n_SubImages
                [~, id] = sort(p_draw(d).corr{sd}, 'descend');
                K = min(sim_draw.nToSave, numel(id));
                for i = 1:K
                    ii = id(i);
                    sim_draw.subimg{sd}(:, i) = p_draw(d).subimg{sd}(:, ii);
                    sim_draw.radius{sd}(i)    = p_draw(d).radius{sd}(ii);
                    sim_draw.x{sd}(i)         = p_draw(d).x{sd}(ii);
                    sim_draw.y{sd}(i)         = p_draw(d).y{sd}(ii);
                    sim_draw.corr{sd}(i)      = p_draw(d).corr{sd}(ii);
                    sim_draw.ssim{sd}(i)      = p_draw(d).ssim{sd}(ii);
                end
                sim_draw.subID{sd} = sd;
                sim_draw.size{sd}  = p_draw(d).size{sd};
            end

            outDir = fullfile(vbl.datadir, vbl.dirList(d).name, 'models');
            if ~exist(outDir, 'dir'), mkdir(outDir); end

            outName = sprintf('%s%s%s.mat', vbl.dirList(d).name, vbl.fileidstr, iff(vbl.debugflag, '_debug', ''));
            save(fullfile(outDir, outName), 'sim_draw');
        end
    end

    function simulate_save_rand_drawings(rep, c, v, trl, tp, vbl)
        % Generate matched "random" phosphenes by reusing the radius and
        % (x, y) coordinates of the saved real singletons but with new
        % cortical map sampling.

        fprintf('Random rep %d / %d\n', rep, vbl.n_Reps);

        % Resample cortical maps per rep (randomness source)
        [c2, v2] = p2p_c.generate_corticalmap(c, v);

        for d = 1:vbl.n_Drawings
            drawDir  = fullfile(vbl.datadir, vbl.dirList(d).name);
            modelsDir = fullfile(drawDir, 'models');
            randDir   = fullfile(drawDir, 'random_models');
            if ~exist(randDir, 'dir'), mkdir(randDir); end

            % Skip if already saved
            tag = iff(vbl.debugflag, '_debug', '');
            outFile = fullfile(randDir, sprintf('%s%s_%d%s.mat', vbl.dirList(d).name, vbl.fileidstr, rep, tag));
            if exist(outFile, 'file'), continue; end

            simFile = fullfile(modelsDir, sprintf('%s%s%s.mat', vbl.dirList(d).name, vbl.fileidstr, tag));
            if ~isfile(simFile)
                warning('Missing sim_draw file for drawing %s; skipping.', vbl.dirList(d).name);
                continue;
            end

            S = load(simFile, 'sim_draw');
            sim_draw = S.sim_draw;

            clear rand_draw
            for sd = 1:sim_draw.n_SubImages
                K = sim_draw.nToSave;
                rand_draw.subimg{sd} = NaN(numel(sim_draw.patient_img{sd}), K);
                rand_draw.corr{sd}   = NaN(K,1);
                rand_draw.ssim{sd}   = NaN(K,1);

                for i = 1:K
                    c2e = c2; v2e = v2;
                    c2e.e.radius = sim_draw.radius{sd}(i);
                    v2e.e.x      = sim_draw.x{sd}(i);
                    v2e.e.y      = sim_draw.y{sd}(i);

                    c2e = p2p_c.define_electrodes(c2e, v2e);
                    c2e = p2p_c.generate_ef(c2e);
                    v2e = p2p_c.generate_corticalelectricalresponse(c2e, v2e);

                    img = ed.generate_phosphene(v2e, tp, trl, vbl);

                    % Align to the corresponding subimage (sd)
                    target = sim_draw.patient_img{sd};
                    [s, r] = findScaleRotationNGC_if(single(img), single(target));
                    [tform, ~] = resolveSimilarityRotationAmbiguityNGC_if(single(img), single(target), s, r);
                    img_aligned = imwarp(img, tform, "OutputView", sim_draw.ref2d);

                    vImg = double(img_aligned(:));
                    vImg = ed.norm255(vImg, vbl.rectify);
                    rand_draw.subimg{sd}(:, i) = vImg(:);

                    rand_draw.corr{sd}(i) = corr(double(target(:)), double(img_aligned(:)));
                    rand_draw.ssim{sd}(i) = ssim(reshape(double(img_aligned), size(target)), double(target), 'Exponents', [0 0 1]);
                end
            end

            saved_ok = false;
            while ~saved_ok
                try
                    save(outFile, 'rand_draw');
                    saved_ok = true;
                catch ME
                    warning('Retry saving random models');
                    pause(0.2);
                end
            end
        end
    end

    % ---------------------------
    % Combining functions
    % ---------------------------
    function combine_sim_draws(vbl)
        warning('off','MATLAB:rankDeficientMatrix');

        for d = 1:vbl.n_Drawings
            drawDir  = fullfile(vbl.datadir, vbl.dirList(d).name);
            modelsDir = fullfile(drawDir, 'models');
            combosDir = fullfile(drawDir, 'combos');
            if ~exist(combosDir,'dir'), mkdir(combosDir); end

            fprintf('Combining real singletons: %d / %d\n', d, vbl.n_Drawings);

            tag = iff(vbl.debugflag, '_debug', '');
            inFile  = fullfile(modelsDir, sprintf('%s%s%s.mat', vbl.dirList(d).name, vbl.fileidstr, tag));
            S = load(inFile, 'sim_draw'); sim_draw = S.sim_draw;

            % Flatten all subimages into columns
            tmp_model = struct('subimg', [], 'corr', [], 'x', [], 'y', [], 'radius', [], 'subID', []);
            ct = 1;
            for sd = 1:sim_draw.n_SubImages
                for ii = 1:sim_draw.nToSave
                    tmp_model.subimg(:, ct) = double(sim_draw.subimg{sd}(:, ii));
                    tmp_model.corr(ct)      = sim_draw.corr{sd}(ii);
                    tmp_model.x(ct)         = sim_draw.x{sd}(ii);
                    tmp_model.y(ct)         = sim_draw.y{sd}(ii);
                    tmp_model.radius(ct)    = sim_draw.radius{sd}(ii);
                    tmp_model.subID(ct)     = sim_draw.subID{sd};
                    ct = ct + 1;
                end
            end

            % All binary combos up to n_DrawCombined
            nCols = size(tmp_model.subimg, 2);
            combos = dec2bin(1:(2^nCols - 1)) == '1';
            combos = combos(sum(combos, 2) <= vbl.n_DrawCombined, :);

            nC = size(combos, 1);
            nimg    = NaN(nC, 1);
            corr_v  = NaN(nC, 1);
            ssim_v  = NaN(nC, 1);
            cmbx    = false(nC, nCols);
            subimg  = NaN(nC, size(tmp_model.subimg, 1));
            target  = double(sim_draw.patient_img{1}(:)); % use first subimage as target for combining

            for ci = 1:nC
                mask = combos(ci, :);
                X    = double(tmp_model.subimg(:, mask)) / 255;
                B    = X \ (target / 255);
                recon = X * B;

                nimg(ci)   = sum(mask);
                corr_v(ci) = corr(recon, target);
                ssim_v(ci) = ssim(reshape(recon, size(sim_draw.patient_img{1})), double(sim_draw.patient_img{1}), 'Exponents', [0 0 1]);
                cmbx(ci, :) = mask;
                subimg(ci, :) = recon;
            end

            [~, id] = sort(corr_v, 'descend');
            combo.nimg     = nimg(id);
            combo.corr_val = corr_v(id);
            combo.ssim_val = ssim_v(id);
            combo.cmbx     = cmbx(id, :);
            combo.subimg   = subimg(id, :);

            outFile = fullfile(combosDir, sprintf('%s%s%s.mat', vbl.dirList(d).name, vbl.fileidstr, tag));
            % -v7.3 if very large
            save(outFile, 'combo', '-v7.3');
        end
    end

    function combine_random_models(vbl)
        warning('off','MATLAB:rankDeficientMatrix');

        for d = 1:vbl.n_Drawings
            drawName = vbl.dirList(d).name;
            drawDir  = fullfile(vbl.datadir, drawName);
            combosDir = fullfile(drawDir, 'combos');
            modelsDir = fullfile(drawDir, 'models');

            fprintf('Combining random models: %d / %d\n', d, vbl.n_Drawings);

            tag = iff(vbl.debugflag, '_debug', '');
            comboFile = fullfile(combosDir, sprintf('%s%s%s.mat', drawName, vbl.fileidstr, tag));
            simFile   = fullfile(modelsDir,  sprintf('%s%s%s.mat', drawName, vbl.fileidstr, tag));

            C = load(comboFile, 'combo'); combo = C.combo;
            S = load(simFile,   'sim_draw'); sim_draw = S.sim_draw;

            % Flatten real subimages as columns
            tmp_model = [];
            ct = 1;
            for sd = 1:sim_draw.n_SubImages
                for ii = 1:sim_draw.nToSave
                    tmp_model.subimg(:, ct) = double(sim_draw.subimg{sd}(:, ii));
                    ct = ct + 1;
                end
            end

            nComb = size(combo.cmbx, 1);
            R = vbl.n_Reps;
            rand_combo.nimg     = NaN(R, nComb);
            rand_combo.corr_val = NaN(R, nComb);
            rand_combo.ssim_val = NaN(R, nComb);
            rand_combo.id       = cell(R, 1);
            rand_combo.perms    = cell(R, 1);
            rand_combo.cmbx     = combo.cmbx;

            for rep = 1:R
                tagR = iff(vbl.debugflag, sprintf('_%d_debug', rep), sprintf('_%d', rep));
                randFile = fullfile(drawDir, 'random_models', sprintf('%s%s%s.mat', drawName, vbl.fileidstr, tagR));

                try
                    D = load(randFile, 'rand_draw');
                    rand_draw = D.rand_draw;
                catch
                    warning('Missing random file for rep %d in %s', rep, drawName);
                    continue;
                end

                % Build randomized order per subimage and flattened matrix
                randimg_local = [];
                perms_rep = cell(1, numel(rand_draw.subimg));
                ct = 1;
                for sd = 1:numel(rand_draw.subimg)
                    nK = numel(rand_draw.corr{sd});
                    shuf_ind = randperm(nK);
                    perms_rep{sd} = shuf_ind;
                    for ii = 1:nK
                        randimg_local(:, ct) = double(rand_draw.subimg{sd}(:, shuf_ind(ii)));
                        ct = ct + 1;
                    end
                end

                target = double(sim_draw.patient_img{1}(:));
                local_corr = NaN(1, nComb);
                local_ssim = NaN(1, nComb);
                local_nimg = NaN(1, nComb);

                for ci = 1:nComb
                    mask = logical(combo.cmbx(ci, :));
                    if ~any(mask)
                        continue;
                    end
                    X  = double(tmp_model.subimg(:, mask) / 255);
                    Xr = double(randimg_local(:, mask) / 255);

                    % Overwrite first column with random (consistent with original)
                    X(:,1) = Xr(:,1);

                    B = X \ target;
                    recon = X * B;

                    local_nimg(ci) = sum(mask);
                    local_corr(ci) = corr(recon, target);
                    local_ssim(ci) = ssim(reshape(recon, size(sim_draw.patient_img{1})), double(sim_draw.patient_img{1}), 'Exponents', [0 0 1]);
                end

                [~, order] = sort(local_corr, 'descend');

                rand_combo.nimg(rep, :)     = local_nimg(order);
                rand_combo.corr_val(rep, :) = local_corr(order);
                rand_combo.ssim_val(rep, :) = local_ssim(order);
                rand_combo.id{rep}          = order;
                rand_combo.perms{rep}       = perms_rep;
            end

            outFile = fullfile(combosDir, sprintf('%s%s%s_rand.mat', drawName, vbl.fileidstr, iff(vbl.debugflag, '_debug', '')));
            save(outFile, 'rand_combo');
        end
    end

    % ---------------------------
    % Visualization
    % ---------------------------
    function visualize_singletons(vbl, d)
        tag = iff(vbl.debugflag, '_debug', '');
        inFile = fullfile(vbl.datadir, vbl.dirList(d).name, 'models', sprintf('%s%s%s.mat', vbl.dirList(d).name, vbl.fileidstr, tag));
        if ~isfile(inFile)
            error('Missing file: %s', inFile);
        end
        S = load(inFile, 'sim_draw'); sim_draw = S.sim_draw;

        figure(1); clf; set(gcf,'Name', vbl.dirList(d).name);
        for i = 1:numel(sim_draw.subimg)
            subplot(2,2,i);
            imagesc(sim_draw.patient_img{i});
            colormap(gray); axis off;
        end

        for s = 1:numel(sim_draw.subimg)
            figure(s+1); clf; set(gcf,'Name', vbl.dirList(d).name);
            K = min(6, size(sim_draw.subimg{s}, 2));
            for i = 1:K
                subplot(3,2,i);
                imagesc(reshape(sim_draw.subimg{s}(:, i), sim_draw.size{s}));
                title(sprintf('corr=%.3f', sim_draw.corr{s}(i)));
                colormap(gray); axis off;
            end
        end
    end

    function visualize_combinations(vbl, d)
        tag = iff(vbl.debugflag, '_debug', '');
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

        for ni = 1:3
            figure(ni+1); clf; set(gcf,'Name', sprintf('Best %d', ni));
            idx = find(combo.nimg == ni);
            K = min(numel(idx), 6);
            for i = 1:K
                subplot(2,3,i);
                imagesc(reshape(combo.subimg(idx(i), :), sim_draw.size{1}));
                title(sprintf('corr=%.3f | idx=%s', combo.corr_val(idx(i)), mat2str(find(combo.cmbx(idx(i), :)))));
                colormap(gray); axis off;
            end
        end
    end

    function show_best_combos(vbl, d)
        % Display target, best real, and best random reconstructions for nimg=1..3
        draw_name = vbl.dirList(d).name;
        tag = iff(vbl.debugflag, '_debug', '');

        combo_file = fullfile(vbl.datadir, draw_name, 'combos', [draw_name, vbl.fileidstr, tag, '.mat']);
        rand_file  = fullfile(vbl.datadir, draw_name, 'combos', [draw_name, vbl.fileidstr, '_rand', tag, '.mat']);
        sim_file   = fullfile(vbl.datadir, draw_name, 'models', [draw_name, vbl.fileidstr, tag, '.mat']);

        C = load(combo_file, 'combo'); combo = C.combo;
        R = load(rand_file, 'rand_combo'); rand_combo = R.rand_combo;
        S = load(sim_file, 'sim_draw'); sim_draw = S.sim_draw;

        target = double(sim_draw.patient_img{1});
        target_vec = target(:);

        tmp_model = [];
        ct = 1;
        for sd = 1:sim_draw.n_SubImages
            for ii = 1:sim_draw.nToSave
                tmp_model.subimg(:, ct) = double(sim_draw.subimg{sd}(:, ii));
                ct = ct + 1;
            end
        end

        figure('Name', sprintf('Drawing %d: %s', d, draw_name), 'Color','w');

        for nimg = 1:3
            recon_real = nan(size(target)); real_corr_val = NaN; real_corr_recomputed = NaN;

            mask_real = (combo.nimg == nimg);
            if any(mask_real)
                [~, best_idx_rel] = max(combo.corr_val(mask_real));
                real_idxs = find(mask_real);
                real_idx = real_idxs(best_idx_rel);
                mask = logical(combo.cmbx(real_idx, :));
                X = double(tmp_model.subimg(:, mask) / 255);
                B = X \ target_vec;
                recon_real = reshape(X * B, size(target));
                real_corr_val = combo.corr_val(real_idx);
                real_corr_recomputed = corr(recon_real(:), target_vec);
            end

            % Best random over reps
            best_corr = -Inf; best_rep = NaN; best_pos = NaN;
            for rep = 1:vbl.n_Reps
                this_nimg = rand_combo.nimg(rep, :);
                this_corr = rand_combo.corr_val(rep, :);
                mask = (this_nimg == nimg);
                if any(mask)
                    [cmax, relpos] = max(this_corr(mask));
                    if cmax > best_corr
                        best_corr = cmax;
                        best_rep = rep;
                        idxs = find(mask);
                        best_pos = idxs(relpos);
                    end
                end
            end

            recon_rand = nan(size(target)); rand_corr_val = NaN; rand_corr_recomputed = NaN;
            if ~isnan(best_rep)
                rand_corr_val = rand_combo.corr_val(best_rep, best_pos);
                rand_c = rand_combo.id{best_rep}(best_pos);
                mask = logical(combo.cmbx(rand_c, :));

                tagR = iff(vbl.debugflag, sprintf('_%d_debug', best_rep), sprintf('_%d', best_rep));
                rand_draw_file = fullfile(vbl.datadir, draw_name, 'random_models', [draw_name, vbl.fileidstr, tagR, '.mat']);
                D = load(rand_draw_file, 'rand_draw');
                rand_draw = D.rand_draw;

                randimg_local = [];
                ct = 1;
                for sd = 1:numel(rand_draw.subimg)
                    shuf_ind = rand_combo.perms{best_rep}{sd};
                    for ii = 1:numel(shuf_ind)
                        randimg_local(:, ct) = double(rand_draw.subimg{sd}(:, shuf_ind(ii)));
                        ct = ct + 1;
                    end
                end

                X  = double(tmp_model.subimg(:, mask) / 255);
                Xr = double(randimg_local(:, mask) / 255);
                if any(mask)
                    X(:,1) = Xr(:,1);
                end
                Br = X \ target_vec;
                recon_rand = reshape(X * Br, size(target));
                rand_corr_recomputed = corr(recon_rand(:), target_vec);
            end

            % Plot row nimg
            subplot(3,3,(nimg-1)*3 + 1);
            imagesc(target); axis image off; colormap gray;
            if nimg == 1, title('Target'); end
            ylabel(sprintf('nimg = %d', nimg));

            subplot(3,3,(nimg-1)*3 + 2);
            imagesc(recon_real); axis image off; colormap gray;
            title(sprintf('Real (saved=%.3f, rec=%.3f)', real_corr_val, real_corr_recomputed));

            subplot(3,3,(nimg-1)*3 + 3);
            imagesc(recon_rand); axis image off; colormap gray;
            title(sprintf('Rand (saved=%.3f, rec=%.3f)', rand_corr_val, rand_corr_recomputed));
        end

        outdir = fullfile(vbl.datadir, 'figures');
        if ~exist(outdir, 'dir'), mkdir(outdir); end
        fig_name = sprintf('Drawing_%02d_%s', d, draw_name);
        exportgraphics(gcf, fullfile(outdir, [fig_name, '.pdf']), 'ContentType', 'vector');
    end

    function plot_corr_histograms(vbl, d)
        draw_name = vbl.dirList(d).name;
        tag = iff(vbl.debugflag, '_debug', '');
        combo_file = fullfile(vbl.datadir, draw_name, 'combos', [draw_name, vbl.fileidstr, tag, '.mat']);
        rand_file  = fullfile(vbl.datadir, draw_name, 'combos', [draw_name, vbl.fileidstr, '_rand', tag, '.mat']);

        C = load(combo_file, 'combo'); combo = C.combo;
        R = load(rand_file,  'rand_combo'); rand_combo = R.rand_combo;

        figure('Name', sprintf('Correlation Histograms: %s', draw_name), 'Color','w');
        for nimg = 1:3
            mask_real = (combo.nimg == nimg);
            best_real_corr = iff(any(mask_real), max(combo.corr_val(mask_real)), NaN);

            best_rand_corrs = NaN(1, vbl.n_Reps);
            for rep = 1:vbl.n_Reps
                this_nimg = rand_combo.nimg(rep, :);
                this_corr = rand_combo.corr_val(rep, :);
                mask = (this_nimg == nimg);
                if any(mask)
                    best_rand_corrs(rep) = max(this_corr(mask));
                end
            end

            subplot(1,3,nimg);
            histogram(best_rand_corrs, 'FaceColor', [0.3 0.3 0.8], 'EdgeColor','k','FaceAlpha',0.6);
            hold on;
            yl = ylim;
            if ~isnan(best_real_corr)
                plot([best_real_corr best_real_corr], yl, 'r-', 'LineWidth', 2);
            end
            ylim(yl); hold off;
            xlabel('Correlation'); ylabel('# Reps');
            title(sprintf('nimg = %d\nReal=%.3f', nimg, best_real_corr));
            legend({'Random bests','Best real'}, 'Location','best');
        end

        outdir = fullfile(vbl.datadir, 'figures');
        if ~exist(outdir, 'dir'), mkdir(outdir); end
        fig_name = sprintf('Drawing_%02d_%s_histogram', d, draw_name);
        exportgraphics(gcf, fullfile(outdir, [fig_name, '.pdf']), 'ContentType', 'vector');
    end

    function compare_real_vs_rand_stats(vbl, d)
        draw_name = vbl.dirList(d).name;
        tag = iff(vbl.debugflag, '_debug', '');
        combo_file = fullfile(vbl.datadir, draw_name, 'combos', [draw_name, vbl.fileidstr, tag, '.mat']);
        rand_file  = fullfile(vbl.datadir, draw_name, 'combos', [draw_name, vbl.fileidstr, '_rand', tag, '.mat']);

        C = load(combo_file, 'combo'); combo = C.combo;
        R = load(rand_file,  'rand_combo'); rand_combo = R.rand_combo;

        for nimg = 1:3
            mask_real = (combo.nimg == nimg);
            best_real_corr = NaN;
            if any(mask_real)
                best_real_corr = max(combo.corr_val(mask_real));
            end

            best_rand_per_rep = NaN(vbl.n_Reps, 1);
            for rep = 1:vbl.n_Reps
                mask = (rand_combo.nimg(rep, :) == nimg);
                if any(mask)
                    best_rand_per_rep(rep) = max(rand_combo.corr_val(rep, mask));
                end
            end
            best_rand_per_rep = best_rand_per_rep(~isnan(best_rand_per_rep));

            top5 = prctile(best_rand_per_rep, 95);
            top1 = prctile(best_rand_per_rep, 99);

            fprintf('Image %s | nimg=%d\n', draw_name, nimg);
            fprintf('Best REAL corr: %.4f\n', best_real_corr);
            fprintf('Random top 5%% threshold: %.4f\n', top5);
            fprintf('Random top 1%% threshold: %.4f\n', top1);

            if best_real_corr >= top1
                fprintf('Real is within TOP 1%% of randoms.\n');
            elseif best_real_corr >= top5
                fprintf('Real is within TOP 5%% of randoms.\n');
            else
                fprintf('Real is BELOW top 5%% of randoms.\n');
            end

            percentile = sum(best_rand_per_rep <= best_real_corr) / numel(best_rand_per_rep) * 100;
            fprintf('Real correlation is at the %.2f percentile\n\n', percentile);
        end
    end

    % ---------------------------
    % Small helpers
    % ---------------------------
    function img = norm01(img)
        img = double(img);
        rng = max(img(:)) - min(img(:));
        if rng < eps
            img = zeros(size(img));
        else
            img = (img - min(img(:))) / rng;
        end
    end

    function vImg = norm255(v, rectify)
        % Normalize vector v to 0..255, optionally rectifying negatives
        v = double(v(:));
        vmax = max(v); vmin = min(v);
        if rectify
            v = max(v - vmin, 0);
            if vmax > 0
                v = v / vmax * 255;
            end
            vImg = uint8(v);
        else
            denom = max(vmax - vmin, eps);
            vImg = int8(128 * (v - vmin) / denom);
        end
    end

    function S = safe_rmfield(S, names)
        if isempty(S), return; end
        if ~isstruct(S), return; end
        keep = isfield(S, names);
        if any(keep)
            S = rmfield(S, names(keep));
        end
    end

    function ensure_dirs(drawDir)
        subdirs = {'models', 'combos', 'random_models', 'drawings'};
        for i = 1:numel(subdirs)
            p = fullfile(drawDir, subdirs{i});
            if ~exist(p, 'dir'), mkdir(p); end
        end
    end

    function clean_dir_safe(p)
        if ~exist(p, 'dir'), return; end
        d = dir(fullfile(p, '*'));
        for i = 1:numel(d)
            if d(i).isdir, continue; end
            try
                delete(fullfile(p, d(i).name));
            catch
            end
        end
    end

end
end

% Simple inline helper
function out = iff(cond, a, b)
if cond, out = a; else, out = b; end
end
