% Simulate_Elche_Drawings_Main.m
% Elche Drawings pipeline
% - Loads patient drawings (and their cropped subimages)
% - Defines cortical/visual model and simulates single-electrode phosphenes
% - Selects "best" singletons per subimage by correlation (and SSIM)
% - Combines singletons into small sets (nimg <= vbl.n_DrawCombined)
% - Generates matched random phosphenes and their combinations
% - Compares real vs random statistics and visualizes results
%
% MATLAB: R2024b tested
% Authors: IF (3/7/2025), comments and functions by ES (7/25/2025, Sep 2025)
% Cleanup, robustness, and documentation pass by ChatGPT (2026)
%
% Requirements on path:
% - ed class (this repo)
% - imregcode helpers: findScaleRotationNGC_if, resolveSimilarityRotationAmbiguityNGC_if
% - p2p_c toolbox: define_* and generate_* functions used throughout

clear; close all; clc;

% ---------------------------
% Setup and configuration
% ---------------------------
% Choose a location token:
%   'auto'  -> infer from computer arch and use current folder as datadir
%   'Z'     -> use your mounted Z: path
%   custom  -> pass a full path if desired
vbl = ed.setup(computer);        % or 'auto' for current folder


% Optional: seed RNG for reproducibility of electrode positions and sampling
% rng(42, 'twister');

% ---------------------------
% Load/crop patient drawings
% ---------------------------
% If you haven’t created cropped versions yet, run once:
% ed.crop_drawings(vbl);

p_draw = ed.load_patient_drawings(vbl);

% ---------------------------
% Define cortical/visual model
% ---------------------------
[c, v, trl, tp] = ed.define_cortical_model(vbl);

% ---------------------------
% Generate and save "best" singletons
% ---------------------------
DO_SIMULATE_AND_COMBINE = true;
if DO_SIMULATE_AND_COMBINE
    p_draw = ed.simulate_drawings(c, v, trl, tp, p_draw, vbl);
    ed.save_simulated_drawings(p_draw, vbl);

    % Visualize singletons for a drawing index (example: 1)
    % ed.visualize_singletons(vbl, 1);

    % Combine singletons into small sets; save combos
    ed.combine_sim_draws(vbl);

    % Optional: visualize combinations for each drawing
    for d = 1:numel(vbl.dirList)
        ed.visualize_combinations(vbl, d);
    end
end

% ---------------------------
% Create the random versions
% ---------------------------
% Remove heavy fields before passing copies to workers
c = ed.safe_rmfield(c, {'e','x','y','X','Y','v','cropPix'});
v = ed.safe_rmfield(v, {'e','x','y','X','Y'});

% Ensure parallel pool up once
pool = gcp('nocreate');
if isempty(pool)
    % Adjust based on your machine
    parpool('threads'); % or parpool(10) for process-based pool
end

% Total reps; each rep saved per drawing
totalReps = vbl.n_Reps;
parfor rep = 1:totalReps
    ed.simulate_save_rand_drawings(rep, c, v, trl, tp, vbl);
end

% Optionally shut down pool after use
% delete(gcp('nocreate'));

disp('Done generating random images');

% ---------------------------
% Combine random models
% ---------------------------
ed.combine_random_models(vbl);

% ---------------------------
% Analysis: real vs random
% ---------------------------
for d = 1:numel(vbl.dirList)
    % ed.show_best_combos(vbl, d);       % target vs best real vs best rand
    % ed.plot_corr_histograms(vbl, d);   % hist of best rand vs best real
    ed.compare_real_vs_rand_stats(vbl, d); % percentiles
end

% Example of targeted plots for last drawing explicitly
lastIdx = numel(vbl.dirList);
ed.plot_corr_histograms(vbl, lastIdx);
ed.compare_real_vs_rand_stats(vbl, lastIdx);
