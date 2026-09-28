
clear; close all;

% Set up path
if strcmp(computer, 'GLNXA64')
    datadir = fullfile(pwd);
else
    datadir = fullfile('Z:', 'IoneFine', 'MyProjects', 'p2p-cortical', 'ElcheDrawings');
end
cd(datadir)
dirList = dir([datadir, '/*img*']);
load(fullfile(datadir, 'patient_draw.mat'));
n_Drawings = length(patient_draw);

% Initialize combined storage
all_real_2 = [];
all_real_3 = [];
all_rand_2 = [];
all_rand_3 = [];

% Initialize per-drawing metrics
metrics = struct([]);

for d = 1:n_Drawings
    % Minimal progress print
    fprintf('Processing drawing %d / %d\n', d, n_Drawings);

    % Load real and random results
    load(fullfile(datadir, dirList(d).name, 'best_combinations.mat'), 'combo');
    load(fullfile(datadir, dirList(d).name, 'rand_combinations.mat'), 'r_combo');

    % Extract real z-values
    real_2 = combo.zval(combo.nimg == 2);
    real_3 = combo.zval(combo.nimg == 3);

    % Extract random z-values
    rand_2 = [];
    rand_3 = [];
    for rep = 1:length(r_combo)
        rand_2 = [rand_2; r_combo(rep).zval(r_combo(rep).nimg == 2)];
        rand_3 = [rand_3; r_combo(rep).zval(r_combo(rep).nimg == 3)];
    end

    % Clean vectors: remove NaN/Inf, force column
    real_2 = real_2(:); rand_2 = rand_2(:);
    real_3 = real_3(:); rand_3 = rand_3(:);
    real_2 = real_2(isfinite(real_2)); rand_2 = rand_2(isfinite(rand_2));
    real_3 = real_3(isfinite(real_3)); rand_3 = rand_3(isfinite(rand_3));

    % Append to full group distributions
    all_real_2 = [all_real_2; real_2];
    all_rand_2 = [all_rand_2; rand_2];
    all_real_3 = [all_real_3; real_3];
    all_rand_3 = [all_rand_3; rand_3];

    % Per-drawing summary
    metrics(d).drawing_id = d;
    metrics(d).real_2_mean = mean(real_2);
    metrics(d).rand_2_mean = mean(rand_2);
    metrics(d).real_3_mean = mean(real_3);
    metrics(d).rand_3_mean = mean(rand_3);

    % Statistical tests for 2-image combos
    if numel(real_2) >= 2 && numel(rand_2) >= 2
        [~, metrics(d).ttest_p_2] = ttest2(real_2, rand_2);
        metrics(d).ranksum_p_2 = ranksum(real_2, rand_2);
    else
        metrics(d).ttest_p_2 = NaN;
        metrics(d).ranksum_p_2 = NaN;
    end

    % Statistical tests for 3-image combos
    if numel(real_3) >= 2 && numel(rand_3) >= 2
        [~, metrics(d).ttest_p_3] = ttest2(real_3, rand_3);
        metrics(d).ranksum_p_3 = ranksum(real_3, rand_3);
    else
        metrics(d).ttest_p_3 = NaN;
        metrics(d).ranksum_p_3 = NaN;
    end
end

%%  Plot histograms of pooled results
figure;
subplot(2,2,1); histogram(all_real_2, 'FaceAlpha', 0.6); title('Real z (nimg = 2)');
xlabel('z'); ylabel('count');
subplot(2,2,2); histogram(all_rand_2, 'FaceAlpha', 0.6); title('Random z (nimg = 2)');
xlabel('z'); ylabel('count');
subplot(2,2,3); histogram(all_real_3, 'FaceAlpha', 0.6); title('Real z (nimg = 3)');
xlabel('z'); ylabel('count');
subplot(2,2,4); histogram(all_rand_3, 'FaceAlpha', 0.6); title('Random z (nimg = 3)');
xlabel('z'); ylabel('count');

%% Create and display summary table
T = struct2table(metrics);
disp(T)

% plot mean z-values per drawing
figure;
subplot(1,2,1);
bar([T.real_2_mean, T.rand_2_mean]);
legend('Real', 'Random');
xlabel('Drawing'); ylabel('z');
title('Mean z (nimg = 2) per drawing');

subplot(1,2,2);
bar([T.real_3_mean, T.rand_3_mean]);
legend('Real', 'Random');
xlabel('Drawing'); ylabel('z');
title('Mean z (nimg = 3) per drawing');