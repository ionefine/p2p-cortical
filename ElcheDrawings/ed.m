% ed.m
%
% Optimized support functions for the Elche Drawings Project.
% Changes:
% - Single precision for big math (X, target, B, recon)
% - Fast correlation (normalized dot product)
% - Optional coarse pre-screen before alignment
% - No combo.subimg stored to reduce disk/memory
% - Thread-based pool friendly (less copying)
% - Robust I/O and safe field handling
%
% Usage: call as static class methods, e.g., vbl = ed.setup('Z');

classdef ed
    methods(Static)

        % ---------------------------
        % Setup
        % ---------------------------
        function vbl = setup(loc)
            % Setup paths and locate data directory
            here = fileparts(mfilename('fullpath'));
            addpath(here);
            addpath(fullfile(here, '..'));
            addpath(fullfile(here, '..', 'imregcode'));

            if nargin < 1
                loc = 'auto';
            end

            switch loc
                case 'Z'
                    vbl.datadir = fullfile('Z:', 'IoneFine', 'MyProjects', 'p2p-cortical', 'ElcheDrawings');
                case 'auto'
                    vbl.datadir = pwd;
                otherwise
                    if isfolder(loc)
                        vbl.datadir = loc;
                    else
                        vbl.datadir = pwd;
                    end
            end

            % Defaults (can be overridden in main)
            vbl.n_SimDraw      = 24;
            vbl.n_DrawCombined = 3;
            vbl.debugflag      = 0;
            vbl.n_Reps         = ed.iff(vbl.debugflag==1, 4, 100);
            vbl.rectify        = 1;
            vbl.fileidstr      = '_10_2_2025';

            % Ask for clean start unless overridden in main
            vbl.cleanstartflag = 1;
            vbl.cleanstartflag = input('Delete all previous simulations? 1=yes, 0=no ... ');   
            
            vbl.dirList = dir(fullfile(vbl.datadir, '*img*'));
            vbl.n_Drawings = numel(vbl.dirList);
            if vbl.n_Drawings == 0
                error('No drawing folders found in %s (pattern "*img*").', vbl.datadir);
            end
            if vbl.debugflag
                vbl.n_Drawings = 10;
            end
        end

        % ---------------------------
        % Patient drawings I/O
        % ---------------------------
        function p_draw = load_patient_drawings(vbl)
            for d = 1:vbl.n_Drawings
                drawDir = fullfile(vbl.datadir, vbl.dirList(d).name);
                disp(['loading ', drawDir]);

                % Optional: clean previous results
                if vbl.cleanstartflag
                    ed.ensure_dirs(drawDir);
                    ed.clean_dir_safe(fullfile(drawDir, 'combos'));
                    ed.clean_dir_safe(fullfile(drawDir, 'models'));
                    ed.clean_dir_safe(fullfile(drawDir, 'random_models'));
                    ed.clean_dir_safe(drawDir); % remove all files at drawing root (subfolders untouched)
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
                        p_draw(d).ref2d = imref2d(size(p_draw(d).patient_img{sd}));
                        p_draw(d).canonicalSize = size(img);
                    else
                        % Enforce invariant: all subimages must match size(sd==1)
                        if ~isequal(size(img), p_draw(d).canonicalSize)
                            error('load_patient_drawings:subimageSizeMismatch', ...
                                ['Subimage sizes differ within drawing %s (d=%d).\n' ...
                                'sd=1 size: [%d %d]\n' ...
                                'sd=%d size: [%d %d]\n'], ...
                                vbl.dirList(d).name, d, ...
                                p_draw(d).canonicalSize(1), p_draw(d).canonicalSize(2), ...
                                sd, size(img,1), size(img,2));
                        end
                    end

                    % Preallocate per subimage top-K pools
                    K = p_draw(d).nToSave;
                    p_draw(d).sim_img{sd} = zeros([size(img), vbl.n_SimDraw], 'uint8'); % for potential storage
                    p_draw(d).corr{sd}    = -inf(K, 1);
                    p_draw(d).ssim{sd}    = -inf(K, 1);
                    p_draw(d).radius{sd}  = NaN(K, 1);
                    p_draw(d).x{sd}       = NaN(K, 1);
                    p_draw(d).y{sd}       = NaN(K, 1);
                    p_draw(d).subID{sd}   = NaN(K, 1);
                    p_draw(d).subimg{sd}  = NaN(numel(img), K); % vectorized
                end
            end
        end

        function crop_drawings(vbl)
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

                imwrite(ed.norm01(orig_img_tmp), fullfile(drawDir, [drawName, '_crop.png']));
                imwrite(ed.norm01(edit_img_tmp), fullfile(drawDir, [drawName, '_edit_crop_0.png']));
            end
        end

        % ---------------------------
        % Cortical/visual model
        % ---------------------------
        function [c, v, trl, tp] = define_cortical_model(vbl)
            disp('Defining cortical/visual model...');

            tp = p2p_c.define_temporalparameters();
            trl.freq = NaN;
            trl      = p2p_c.define_trial(tp, trl);

            if vbl.debugflag
                v.pixperdeg = 4;  c.pixpermm  = 4;   nElectrodes = 20;
            else
                v.pixperdeg = 10; c.pixpermm  = 10;  nElectrodes = 1000;
            end

            c.cortexHeight = [-30, 30];
            c.cortexLength = [-55, 55];
            c.onoff_ratio  = 0.9;

            v.visfieldHeight = [-30, 30];
            v.visfieldWidth  = [-30, 30];

            v = p2p_c.define_visualmap(v);
            c = p2p_c.define_cortex(c);
            [c, v] = p2p_c.generate_corticalmap(c, v);

            theta = 2*pi*rand(1, nElectrodes);
            r     = exp(4*rand(1, nElectrodes)); r = r*3./max(r(:));
            ecc   = exp((log(15) * rand(1, nElectrodes)));
            c.I_k = 1000;

            [x, y] = pol2cart(theta(:), ecc(:));
            iperm  = randperm(numel(x));
            for i = 1:numel(x)
                v.e(i).x       = x(iperm(i));
                v.e(i).y       = y(iperm(i));
                c.e(i).radius  = r(iperm(i));
            end

            c = p2p_c.define_electrodes(c, v);
            disp(['Defined ', num2str(numel(x)), ' electrodes.']);
        end

        % ---------------------------
        % Simulation of singletons
        % ---------------------------
        function p_draw = simulate_drawings(c_orig, v_orig, trl, tp, p_draw, vbl)
            %SIMULATE_DRAWINGS
            % For each electrode:
            %   1) Generate phosphene (model)
            %   2) For each drawing/subimage: align model to that subimage
            %   3) Keep top-K models per subimage (by correlation), with diversity constraint
            %
            % Key invariant enforced:
            %   size(img_aligned) MUST equal size(target) for the current sd.
            %   If not, throw an error immediately (fail fast).
            %
            % Parallelized across electrodes: each electrode's phosphene
            % generation (p2p_c.generate_corticalelectricalresponse, the
            % dominant cost -- profiling showed ~130s/electrode at
            % production resolution, i.e. tens of hours for a full
            % 1000-electrode run when done serially) plus registration
            % against every drawing/subimage is fully independent of every
            % other electrode, so phase 1 below runs it in a parfor. Phase 2
            % (choosing which pool slot each candidate replaces) reads and
            % mutates p_draw(d).corr{sd}/subimg{sd}, which is genuinely
            % order-dependent (the "quick gate" and diversity check compare
            % against whatever's already been saved), so it stays a plain,
            % strictly-sequential for loop over electrodes in original
            % (1:nElect) order, applied to the phase-1 results -- this
            % reproduces byte-identical output to the old fully-serial
            % version, just with the expensive independent work moved onto
            % multiple workers.

            corrThresh = 0.90; % max allowed similarity among saved models (diversity constraint)
            nElect = numel(c_orig.e);

            % Flatten (drawing, subimage) pairs once; every electrode is
            % registered against the same list of targets.
            targets = struct('img', {}, 'targetSize', {}, 'ref2d', {}, 'd', {}, 'sd', {});
            t = 0;
            for d = 1:vbl.n_Drawings
                for sd = 1:p_draw(d).n_SubImages
                    t = t + 1;
                    targets(t).img = p_draw(d).patient_img{sd};
                    targets(t).targetSize = size(targets(t).img);
                    targets(t).ref2d = p_draw(d).ref2d;
                    targets(t).d = d;
                    targets(t).sd = sd;
                end
            end
            nTargets = numel(targets);

            % -----------------------
            % Phase 1: parallel, per-electrode candidate computation
            % (no shared mutable state -- safe to run in any order)
            % -----------------------
            candOut = cell(nElect, 1);
            parfor eIdx = 1:nElect
                tElec = tic; % local to this iteration only -- no cross-worker state

                c = ed.safe_rmfield(c_orig, {'e'});
                v = ed.safe_rmfield(v_orig, {'e'});
                c.e = c_orig.e(eIdx);
                v.e = v_orig.e(eIdx);

                c = p2p_c.generate_ef(c);
                v = p2p_c.generate_corticalelectricalresponse(c, v);
                img = uint8(ed.generate_phosphene(v, tp, trl, vbl));

                cand = struct('peakcorr', cell(1, nTargets), 'scaled_vec', cell(1, nTargets), ...
                    'newCorr', cell(1, nTargets), 'newSSIM', cell(1, nTargets), 'targetSize', cell(1, nTargets));

                for tt = 1:nTargets
                    target = targets(tt).img;
                    targetSize = targets(tt).targetSize;

                    [s, r] = findScaleRotationNGC_if(single(img), single(target));
                    [tform, peakcorr] = resolveSimilarityRotationAmbiguityNGC_if(single(img), single(target), s, r);
                    img_aligned = imwarp(img, tform, 'OutputView', targets(tt).ref2d, 'FillValues', 0);

                    if ~isequal(size(img_aligned), targetSize)
                        error('simulate_drawings:sizeMismatch', ...
                            'Warp output size != target size (eIdx=%d, t=%d).', eIdx, tt);
                    end

                    pixCount = numel(target);
                    scaled_vec = ed.norm255(double(img_aligned(:)), vbl.rectify); % uint8 (rectify=1) or int8
                    scaled_vec = double(scaled_vec(:)); % store as double for compatibility with existing code

                    if numel(scaled_vec) ~= pixCount
                        error('simulate_drawings:vectorLengthMismatch', ...
                            'Vector length mismatch after warp: expected %d, got %d (t=%d, e=%d).', ...
                            pixCount, numel(scaled_vec), tt, eIdx);
                    end

                    newCorr = ed.fastcorr(single(target(:)), single(img_aligned(:)));
                    newSSIM = ssim( ...
                        reshape(single(img_aligned), targetSize), ...
                        single(target), ...
                        'Exponents', [0 0 1], ...
                        'DynamicRange', 255);

                    cand(tt).peakcorr = peakcorr;
                    cand(tt).scaled_vec = scaled_vec;
                    cand(tt).newCorr = double(newCorr);
                    cand(tt).newSSIM = double(newSSIM);
                    cand(tt).targetSize = targetSize;
                end

                candOut{eIdx} = struct('cand', cand, 'radius', c.e.radius, 'x', v.e.x, 'y', v.e.y);
                fprintf('Electrode %d / %d computed (%.1fs)\n', eIdx, nElect, toc(tElec));
            end

            % -----------------------
            % Phase 2: sequential merge, in original electrode order
            % (this is where order-dependent pool/diversity state lives)
            % -----------------------
            for eIdx = 1:nElect
                fprintf('Electrode %d / %d\n', eIdx, nElect);
                radius = candOut{eIdx}.radius;
                ex = candOut{eIdx}.x;
                ey = candOut{eIdx}.y;
                cand = candOut{eIdx}.cand;

                for tt = 1:nTargets
                    d = targets(tt).d;
                    sd = targets(tt).sd;
                    targetSize = cand(tt).targetSize;
                    peakcorr = cand(tt).peakcorr;

                    % Quick gate: only consider if it can beat the worst saved corr
                    % (Use peakcorr from the registration as a cheap gate.)
                    if peakcorr <= min(p_draw(d).corr{sd})
                        continue;
                    end

                    % Ensure pool matrix height matches target pixels (otherwise upstream is already corrupted)
                    pixCount = prod(targetSize);
                    if size(p_draw(d).subimg{sd}, 1) ~= pixCount
                        error('simulate_drawings:poolSizeMismatch', ...
                            ['Pool pixel dimension mismatch BEFORE write.\n' ...
                            'Drawing (d=%d), sd=%d\n' ...
                            'Expected pool rows=%d, got %d.\n' ...
                            'This indicates earlier writes used the wrong OutputView.'], ...
                            d, sd, pixCount, size(p_draw(d).subimg{sd}, 1));
                    end

                    scaled_vec = cand(tt).scaled_vec;
                    newCorr = cand(tt).newCorr;
                    newSSIM = cand(tt).newSSIM;

                    % Diversity check vs existing saved candidates
                    pool = p_draw(d).subimg{sd}; % [pixCount x K]
                    K = size(pool, 2);

                    if all(~isfinite(p_draw(d).corr{sd})) || all(isnan(pool(:)))
                        maxPoolCorr = 0;
                        mostSimilarIdx = 1;
                    else
                        pc = -inf(1, K, 'single');
                        for k = 1:K
                            if all(isfinite(pool(:,k)))
                                pc(k) = ed.fastcorr(single(pool(:,k)), single(scaled_vec));
                            end
                        end
                        [maxPoolCorr, mostSimilarIdx] = max(pc);
                        if ~isfinite(maxPoolCorr)
                            maxPoolCorr = 0;
                            mostSimilarIdx = 1;
                        end
                    end

                    % Choose which slot to replace
                    if maxPoolCorr > corrThresh
                        ridx = mostSimilarIdx;           % too similar -> replace most similar
                    else
                        [~, ridx] = min(p_draw(d).corr{sd}); % otherwise replace worst corr
                    end

                    % Save if better than what it's replacing
                    if newCorr > p_draw(d).corr{sd}(ridx)
                        p_draw(d).subimg{sd}(:, ridx) = scaled_vec;
                        p_draw(d).radius{sd}(ridx)   = radius;
                        p_draw(d).x{sd}(ridx)        = ex;
                        p_draw(d).y{sd}(ridx)        = ey;
                        p_draw(d).corr{sd}(ridx)     = newCorr;
                        p_draw(d).ssim{sd}(ridx)     = newSSIM;
                        p_draw(d).subID{sd}(ridx)    = sd;
                        p_draw(d).size{sd}           = targetSize;
                    end
                end
            end
        end


        function img = generate_phosphene(v, tp, trl, vbl)
            trl = p2p_c.generate_phosphene(v, tp, trl);
            img = max(trl.max_phosphene, [], 3);
            img = crop_img(img, 20);
            m = max(abs(img(:)));
            if m < eps, m = 1; end
            img = img ./ m;

            if vbl.rectify
                img(img < 0) = 0;
                img = uint8(255 * img);
            else
                img = uint8(127 * (img + 0.5));
            end
        end

        function save_simulated_drawings(p_draw, vbl)
            for d = 1:vbl.n_Drawings
                
                sim_draw = struct();  % or: clear sim_draw
                
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

                tag = ed.iff(vbl.debugflag, '_debug', '');
                outName = sprintf('%s%s%s.mat', vbl.dirList(d).name, vbl.fileidstr, tag);
                save(fullfile(outDir, outName), 'sim_draw');
            end
        end

        % ---------------------------
        % Random phosphenes
        % ---------------------------
        function simulate_save_rand_drawings(rep, c, v, trl, tp, vbl)
            %SIMULATE_SAVE_RAND_DRAWINGS
            % For each drawing:
            %   - Load the saved real sim_draw (to get target + electrode params per slot)
            %   - Resample cortical maps for this rep
            %   - For each subimage and each saved slot: generate a "random" phosphene
            %     using same (radius, x, y) but new cortical maps, align to *that subimage*
            %   - Save rand_draw for this rep
            %
            % Key invariant enforced:
            %   size(img_aligned) MUST equal size(target) (per sd).
            %   If not, throw an error immediately.

            fprintf('Random rep %d / %d\n', rep, vbl.n_Reps);

            % Resample cortical maps per rep (the randomness source)
            [c2, v2] = p2p_c.generate_corticalmap(c, v);

            tag = ed.iff(vbl.debugflag, '_debug', '');

            for d = 1:vbl.n_Drawings
                drawName  = vbl.dirList(d).name;
                
                drawDir   = fullfile(vbl.datadir, drawName);
                modelsDir = fullfile(drawDir, 'models');
                randDir   = fullfile(drawDir, 'random_models');
                if ~exist(randDir, 'dir'), mkdir(randDir); end

                outFile = fullfile(randDir, sprintf('%s%s_%d%s.mat', drawName, vbl.fileidstr, rep, tag));
                if exist(outFile, 'file')
                    continue; % already saved for this drawing+rep
                end

                simFile = fullfile(modelsDir, sprintf('%s%s%s.mat', drawName, vbl.fileidstr, tag));
                if ~isfile(simFile)
                    error('simulate_save_rand_drawings:missingSimDraw', ...
                        'Missing sim_draw file: %s', simFile);
                end

                S = load(simFile, 'sim_draw');
                sim_draw = S.sim_draw;

                clear rand_draw

                for sd = 1:sim_draw.n_SubImages
                    disp([drawName, ' subdrawing ', num2str(sd)]);
                    target = sim_draw.patient_img{sd};
                    targetSize = size(target);
                    pixCount = numel(target);
                    K = sim_draw.nToSave;

                    % Preallocate with correct dimensions for this sd
                    rand_draw.subimg{sd} = NaN(pixCount, K);
                    rand_draw.corr{sd}   = NaN(K, 1);
                    rand_draw.ssim{sd}   = NaN(K, 1);

                    % Output canvas MUST be the current sd size
                    outRef = imref2d(targetSize);

                    for i = 1:K
                        % Use the (radius, x, y) from the real saved slot
                        c2e = c2; v2e = v2;
                        c2e.e.radius = sim_draw.radius{sd}(i);
                        v2e.e.x      = sim_draw.x{sd}(i);
                        v2e.e.y      = sim_draw.y{sd}(i);

                        c2e = p2p_c.define_electrodes(c2e, v2e);
                        c2e = p2p_c.generate_ef(c2e);
                        v2e = p2p_c.generate_corticalelectricalresponse(c2e, v2e);

                        img = uint8(ed.generate_phosphene(v2e, tp, trl, vbl));

                        % Align to THIS subimage sd (not sd=1)
                        [s, r] = findScaleRotationNGC_if(single(img), single(target));

                        % Preserve your prior behavior if desired (optional):
                        % If your old code intentionally skipped rotation, uncomment:
                         r = 0;

                        [tform, ~] = resolveSimilarityRotationAmbiguityNGC_if(single(img), single(target), s, r);

                        img_aligned = imwarp(img, tform, 'OutputView', outRef, 'FillValues', 0);

                        % HARD SIZE CHECK
                        if ~isequal(size(img_aligned), targetSize)
                            error('simulate_save_rand_drawings:sizeMismatch', ...
                                ['Warped random size mismatch.\n' ...
                                'Drawing: %s (d=%d), sd=%d, rep=%d, i=%d\n' ...
                                'Expected target size: [%d %d]\n' ...
                                'Got aligned size:     [%d %d]\n'], ...
                                drawName, d, sd, rep, i, ...
                                targetSize(1), targetSize(2), size(img_aligned,1), size(img_aligned,2));
                        end

                        % Vectorize and normalize to 0..255-ish (consistent with real pool)
                        scaled_vec = ed.norm255(double(img_aligned(:)), vbl.rectify);
                        scaled_vec = double(scaled_vec(:));

                        if numel(scaled_vec) ~= pixCount
                            error('simulate_save_rand_drawings:vectorLengthMismatch', ...
                                'Vector length mismatch after warp: expected %d, got %d (d=%d, sd=%d, rep=%d, i=%d).', ...
                                pixCount, numel(scaled_vec), d, sd, rep, i);
                        end

                        rand_draw.subimg{sd}(:, i) = scaled_vec;

                        % Metrics
                        rand_draw.corr{sd}(i) = double(ed.fastcorr(single(target(:)), single(img_aligned(:))));
                        rand_draw.ssim{sd}(i) = double(ssim( ...
                            reshape(single(img_aligned), targetSize), ...
                            single(target), ...
                            'Exponents', [0 0 1], ...
                            'DynamicRange', 255));
                    end
                end

                % Save (retry loop to be robust on networked filesystems)
                saved_ok = false;
                while ~saved_ok
                    try
                        save(outFile, 'rand_draw');
                        saved_ok = true;
                    catch ME
                        warning('simulate_save_rand_drawings:saveRetry', ...
                            'Trouble saving %s (%s). Retrying...', outFile, ME.message);
                        pause(0.2);
                    end
                end
            end
        end


        % ---------------------------
        % Combine real and random
        % ---------------------------
        function combine_sim_draws(vbl)
            warning('off','MATLAB:rankDeficientMatrix');

            for d = 1:vbl.n_Drawings
                drawDir   = fullfile(vbl.datadir, vbl.dirList(d).name);
                modelsDir = fullfile(drawDir, 'models');
                combosDir = fullfile(drawDir, 'combos');
                if ~exist(combosDir,'dir'), mkdir(combosDir); end

                fprintf('Combining real singletons: %d / %d\n', d, vbl.n_Drawings);

                tag = ed.iff(vbl.debugflag, '_debug', '');
                inFile  = fullfile(modelsDir, sprintf('%s%s%s.mat', vbl.dirList(d).name, vbl.fileidstr, tag));
                S = load(inFile, 'sim_draw'); sim_draw = S.sim_draw;

                % Flatten all subimages into columns
                tmp_model = struct('subimg', []);
                ct = 1;
                for sd = 1:sim_draw.n_SubImages
                    for ii = 1:sim_draw.nToSave
                        tmp_model.subimg(:, ct) = double(sim_draw.subimg{sd}(:, ii));
                        ct = ct + 1;
                    end
                end

                % Enumerate all combinations up to vbl.n_DrawCombined
                nCols  = size(tmp_model.subimg, 2);
                combos = dec2bin(1:(2^nCols - 1)) == '1';
                combos = combos(sum(combos, 2) <= vbl.n_DrawCombined, :);

                nC     = size(combos, 1);
                nimg    = NaN(nC, 1, 'single');
                corr_v  = NaN(nC, 1, 'single');
                ssim_v  = NaN(nC, 1, 'single');
                cmbx    = false(nC, nCols);
                targetImg = single(sim_draw.patient_img{1});
                target    = targetImg(:);

                for ci = 1:nC
                    mask  = combos(ci, :);
                    X     = single(tmp_model.subimg(:, mask)) / 255;
                    B     = X \ target;
                    recon = X * B;
                    nimg(ci)   = sum(mask);
                    corr_v(ci) = ed.fastcorr(recon, target);
                    ssim_v(ci) = ssim(reshape(recon, size(targetImg)), targetImg, ...
                        'Exponents', [0 0 1], 'DynamicRange', 255);
                    cmbx(ci, :) = mask;
                end

                [~, id] = sort(corr_v, 'descend');
                combo.nimg     = double(nimg(id));
                combo.corr_val = double(corr_v(id));
                combo.ssim_val = double(ssim_v(id));
                combo.cmbx     = cmbx(id, :);
                % Note: combo.subimg omitted to reduce disk footprint

                outFile = fullfile(combosDir, sprintf('%s%s%s.mat', vbl.dirList(d).name, vbl.fileidstr, tag));
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

                tag = ed.iff(vbl.debugflag, '_debug', '');
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
                rand_combo.nimg     = NaN(R, nComb, 'single');
                rand_combo.corr_val = NaN(R, nComb, 'single');
                rand_combo.ssim_val = NaN(R, nComb, 'single');
                rand_combo.id       = cell(R, 1);
                rand_combo.perms    = cell(R, 1);
                rand_combo.cmbx     = combo.cmbx;

                targetImg = single(sim_draw.patient_img{1});
                target    = targetImg(:);

                for rep = 1:R
                    tagR = ed.iff(vbl.debugflag, sprintf('_%d_debug', rep), sprintf('_%d', rep));
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

                    local_corr = NaN(1, nComb, 'single');
                    local_ssim = NaN(1, nComb, 'single');
                    local_nimg = NaN(1, nComb, 'single');

                    for ci = 1:nComb
                        mask = logical(combo.cmbx(ci, :));
                        if ~any(mask), continue; end
                        X  = single(tmp_model.subimg(:, mask)) / 255;
                        Xr = single(randimg_local(:, mask)) / 255;

                        % Overwrite first column with random
                        X(:,1) = Xr(:,1);

                        B = X \ target;
                        recon = X * B;

                        local_nimg(ci) = sum(mask);
                        local_corr(ci) = ed.fastcorr(recon, target);
                        local_ssim(ci) = ssim(reshape(recon, size(targetImg)), targetImg, ...
                            'Exponents', [0 0 1], 'DynamicRange', 255);
                    end

                    [~, order] = sort(local_corr, 'descend');

                    rand_combo.nimg(rep, :)     = local_nimg(order);
                    rand_combo.corr_val(rep, :) = local_corr(order);
                    rand_combo.ssim_val(rep, :) = local_ssim(order);
                    rand_combo.id{rep}          = order;
                    rand_combo.perms{rep}       = perms_rep;
                end

                outFile = fullfile(combosDir, sprintf('%s%s_rand%s.mat', drawName, vbl.fileidstr, tag));
                save(outFile, 'rand_combo', '-v7.3');
            end
        end

        % ---------------------------
        % Visualization
        % ---------------------------
        function visualize_singletons(vbl, d)
            tag = ed.iff(vbl.debugflag, '_debug', '');
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

            for ni = 1:3
                figure(ni+1); clf; set(gcf,'Name', sprintf('Best %d', ni));
                idx = find(combo.nimg == ni);
                K = min(numel(idx), 6);
                for i = 1:K
                    mask  = logical(combo.cmbx(idx(i), :));
                    X     = single(tmp_model.subimg(:, mask)) / 255;
                    B     = X \ target;
                    recon = reshape(X * B, size(targetImg));
                    subplot(2,3,(i));
                    imagesc(recon); axis image off; colormap gray;
                    title(sprintf('corr=%.3f | idx=%s', combo.corr_val(idx(i)), mat2str(find(mask))));
                end
            end
        end

        function show_best_combos(vbl, d)
            draw_name = vbl.dirList(d).name;
            tag = ed.iff(vbl.debugflag, '_debug', '');

            combo_file = fullfile(vbl.datadir, draw_name, 'combos', [draw_name, vbl.fileidstr, tag, '.mat']);
            rand_file  = fullfile(vbl.datadir, draw_name, 'combos', [draw_name, vbl.fileidstr, '_rand', tag, '.mat']);
            sim_file   = fullfile(vbl.datadir, draw_name, 'models', [draw_name, vbl.fileidstr, tag, '.mat']);

            C = load(combo_file, 'combo'); combo = C.combo;
            R = load(rand_file, 'rand_combo'); rand_combo = R.rand_combo;
            S = load(sim_file, 'sim_draw'); sim_draw = S.sim_draw;

            targetImg = double(sim_draw.patient_img{1});
            target_vec = targetImg(:);

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
                recon_real = nan(size(targetImg)); real_corr_val = NaN; real_corr_recomputed = NaN;

                mask_real = (combo.nimg == nimg);
                if any(mask_real)
                    [~, best_idx_rel] = max(combo.corr_val(mask_real));
                    real_idxs = find(mask_real);
                    real_idx = real_idxs(best_idx_rel);
                    mask = logical(combo.cmbx(real_idx, :));
                    X = single(tmp_model.subimg(:, mask)) / 255;
                    B = X \ single(target_vec);
                    recon_real = reshape(X * B, size(targetImg));
                    real_corr_val = combo.corr_val(real_idx);
                    real_corr_recomputed = ed.fastcorr(single(recon_real(:)), single(target_vec));
                end

                % Best random across reps
                best_corr = -Inf; best_rep = NaN; best_pos = NaN;
                for rep = 1:vbl.n_Reps
                    this_nimg = rand_combo.nimg(rep, :);
                    this_corr = rand_combo.corr_val(rep, :);
                    mask = (this_nimg == nimg);
                    if any(mask)
                        [cmax, relpos] = max(this_corr(mask));
                        if cmax > best_corr
                            best_corr = cmax; best_rep = rep;
                            idxs = find(mask);
                            best_pos = idxs(relpos);
                        end
                    end
                end

                recon_rand = nan(size(targetImg)); rand_corr_val = NaN; rand_corr_recomputed = NaN;
                if ~isnan(best_rep)
                    rand_corr_val = rand_combo.corr_val(best_rep, best_pos);
                    rand_c = rand_combo.id{best_rep}(best_pos);
                    mask = logical(combo.cmbx(rand_c, :));

                    tagR = ed.iff(vbl.debugflag, sprintf('_%d_debug', best_rep), sprintf('_%d', best_rep));
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

                    X  = single(tmp_model.subimg(:, mask)) / 255;
                    Xr = single(randimg_local(:, mask)) / 255;
                    if any(mask)
                        X(:,1) = Xr(:,1);
                    end
                    Br = X \ single(target_vec);
                    recon_rand = reshape(X * Br, size(targetImg));
                    rand_corr_recomputed = ed.fastcorr(single(recon_rand(:)), single(target_vec));
                end

                % Plot
                subplot(3,3,(nimg-1)*3 + 1);
                imagesc(targetImg); axis image off; colormap gray;
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
            fig_name = sprintf('Drawing_%s', draw_name);
            exportgraphics(gcf, fullfile(outdir, [fig_name, '.pdf']), 'ContentType', 'vector');
        end

        function plot_corr_histograms(vbl, d)
            draw_name = vbl.dirList(d).name;
            tag = ed.iff(vbl.debugflag, '_debug', '');
            combo_file = fullfile(vbl.datadir, draw_name, 'combos', [draw_name, vbl.fileidstr, tag, '.mat']);
            rand_file  = fullfile(vbl.datadir, draw_name, 'combos', [draw_name, vbl.fileidstr, '_rand', tag, '.mat']);

            C = load(combo_file, 'combo'); combo = C.combo;
            R = load(rand_file,  'rand_combo'); rand_combo = R.rand_combo;

            figure('Name', sprintf('Correlation Histograms: %s', draw_name), 'Color','w');
            for nimg = 1:3
                mask_real = (combo.nimg == nimg);
                best_real_corr = ed.iff(any(mask_real), max(combo.corr_val(mask_real)), NaN);

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
            fig_name = sprintf('Drawing_%s_histogram', draw_name);
            exportgraphics(gcf, fullfile(outdir, [fig_name, '.pdf']), 'ContentType', 'vector');
        end

        function compare_real_vs_rand_stats(vbl, d)
            draw_name = vbl.dirList(d).name;
            tag = ed.iff(vbl.debugflag, '_debug', '');
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
            if isempty(S) || ~isstruct(S), return; end
            names = names(isfield(S, names));
            if ~isempty(names), S = rmfield(S, names); end
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
                try, delete(fullfile(p, d(i).name)); catch, end
            end
        end

        function r = fastcorr(a, b)
            % Fast NCC: correlation via normalized dot product
            a = single(a(:)); b = single(b(:));
            a = a - mean(a); b = b - mean(b);
            na = sqrt(sum(a.*a)); nb = sqrt(sum(b.*b));
            if na == 0 || nb == 0
                r = single(0);
            else
                r = sum(a.*b) / (na*nb);
            end
        end

        function out = iff(cond, a, b)
            if cond, out = a; else, out = b; end
        end

    end
end
