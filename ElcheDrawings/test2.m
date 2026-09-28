clear; close all

topNum = 1; 
topNum2 = 1;

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

% Initialize 
all_real_1 = [];
all_real_2 = [];
all_real_3 = [];
all_rand_1 = [];
all_rand_2 = [];
all_rand_3 = [];

metrics = struct([]);

for d = 1:n_Drawings
    fprintf('Processing drawing %d / %d\n', d, n_Drawings);

    load(fullfile(datadir, dirList(d).name, 'best_combinations.mat'), 'combo');
    load(fullfile(datadir, dirList(d).name, 'rand_combinations.mat'), 'r_combo');

    % Extract z-values for each nimg
    real_1 = combo.zval(combo.nimg == 1);
    real_2 = combo.zval(combo.nimg == 2);
    real_3 = combo.zval(combo.nimg == 3);

    rand_1 = []; rand_2 = []; rand_3 = [];
    for rep = 1:length(r_combo)
        z1 = r_combo(rep).zval(r_combo(rep).nimg ==1);
        %z1 = randperm()sort(z1, 'descend');
        z1_rand = z1(randperm(length(z1), 1));
        rand_1 = [rand_1, z1_rand];
        %rand_1 = [rand_1; z1(1:min(topNum2, numel(z1)))];

        z2 = r_combo(rep).zval(r_combo(rep).nimg ==2);
        z2_rand = z2(randperm(length(z2), 1));
        rand_2 = [rand_2, z2_rand];

        %z2 = sort(z2, 'descend')
        %rand_2 = [rand_2; z2(1:min(topNum2, numel(z2)))];


        z3 = r_combo(rep).zval(r_combo(rep).nimg ==3);
        %z3= sort(z3, 'descend');
        %rand_3 = [rand_3; z3(1:min(topNum2, numel(z3)))];

        z3_rand = z3(randperm(length(z3), 1));
        rand_3 = [rand_3, z3_rand];

        %rand_1 = [rand_1; r_combo(rep).zval(r_combo(rep).nimg == 1)];
        %rand_2 = [rand_2; r_combo(rep).zval(r_combo(rep).nimg == 2)];
        %rand_3 = [rand_3; r_combo(rep).zval(r_combo(rep).nimg == 3)];
    end

    % Clean 
    real_1 = real_1(:); real_2 = real_2(:); real_3 = real_3(:);
    rand_1 = rand_1(:); rand_2 = rand_2(:); rand_3 = rand_3(:);
    real_1 = real_1(isfinite(real_1)); rand_1 = rand_1(isfinite(rand_1));
    real_3 = real_3(isfinite(real_3)); rand_3 = rand_3(isfinite(rand_3));

    real_2 = real_2(isfinite(real_2)); rand_2 = rand_2(isfinite(rand_2));
    % Cap both real and random to topNum highest values
    real_1 = sort(real_1, 'descend'); if numel(real_1) > topNum, real_1 = real_1(1:topNum); end
    real_2 = sort(real_2, 'descend'); if numel(real_2) > topNum, real_2 = real_2(1:topNum); end
    real_3 = sort(real_3, 'descend'); if numel(real_3) > topNum, real_3 = real_3(1:topNum); end
    

    % Save for plotting % new histogram
    all_top1_rand_per_drawing{d} = rand_1;  % already top-1 per rep
    all_top1_real_per_drawing(d) = real_1;  % already top-1 overall

    %rand_1 = sort(rand_1, 'descend'); if numel(rand_1) > topNum, rand_1 = rand_1(1:topNum); end
    %rand_2 = sort(rand_2, 'descend'); if numel(rand_2) > topNum, rand_2 = rand_2(1:topNum); end
    %rand_3 = sort(rand_3, 'descend'); if numel(rand_3) > topNum, rand_3 = rand_3(1:topNum); end

    
    %  Cap randoms to match real counts 
    %N1 = numel(real_1); N2 = numel(real_2); N3 = numel(real_3);
    %rand_1 = sort(rand_1, 'descend'); if numel(rand_1) > N1, rand_1 = rand_1(1:N1); end
    % rand_2 = sort(rand_2, 'descend'); if numel(rand_2) > N2, rand_2 = rand_2(1:N2); end
    %rand_3 = sort(rand_3, 'descend'); if numel(rand_3) > N3, rand_3 = rand_3(1:N3); end

    % Append to group distributions 
    all_real_1 = [all_real_1; real_1];
    all_rand_1 = [all_rand_1; rand_1];
    all_real_2 = [all_real_2; real_2];
    all_rand_2 = [all_rand_2; rand_2];
    all_real_3 = [all_real_3; real_3];
    all_rand_3 = [all_rand_3; rand_3];

    % Per-drawing stats 
    metrics(d).drawing_id = d;
    metrics(d).real_1_mean = mean(real_1);
    metrics(d).rand_1_mean = mean(rand_1);
    metrics(d).real_2_mean = mean(real_2);
    metrics(d).rand_2_mean = mean(rand_2);
    metrics(d).real_3_mean = mean(real_3);
    metrics(d).rand_3_mean = mean(rand_3);

    % Statistical tests
    if numel(real_1) >= 2 && numel(rand_1) >= 2
        [~, metrics(d).ttest_p_1] = ttest2(real_1, rand_1);
        metrics(d).ranksum_p_1 = ranksum(real_1, rand_1);
    else
        metrics(d).ttest_p_1 = NaN;
        metrics(d).ranksum_p_1 = NaN;
    end
    if numel(real_2) >= 2 && numel(rand_2) >= 2
        [~, metrics(d).ttest_p_2] = ttest2(real_2, rand_2);
        metrics(d).ranksum_p_2 = ranksum(real_2, rand_2);
    else
        metrics(d).ttest_p_2 = NaN;
        metrics(d).ranksum_p_2 = NaN;
    end
    if numel(real_3) >= 2 && numel(rand_3) >= 2
        [~, metrics(d).ttest_p_3] = ttest2(real_3, rand_3);
        metrics(d).ranksum_p_3 = ranksum(real_3, rand_3);
    else
        metrics(d).ttest_p_3 = NaN;
        metrics(d).ranksum_p_3 = NaN;
    end
end

% Overall statistics
overall.real_1_mean = mean(all_real_1);
overall.rand_1_mean = mean(all_rand_1);
overall.real_2_mean = mean(all_real_2);
overall.rand_2_mean = mean(all_rand_2);
overall.real_3_mean = mean(all_real_3);
overall.rand_3_mean = mean(all_rand_3);

% t-tests
[~, overall.ttest_p_1] = ttest2(all_real_1, all_rand_1);
[~, overall.ttest_p_2] = ttest2(all_real_2, all_rand_2);
[~, overall.ttest_p_3] = ttest2(all_real_3, all_rand_3);

% rank-sum (non-parametric)
overall.ranksum_p_1 = ranksum(all_real_1, all_rand_1);
overall.ranksum_p_2 = ranksum(all_real_2, all_rand_2);
overall.ranksum_p_3 = ranksum(all_real_3, all_rand_3);

% Display
disp('OVERALL STATS')
fprintf('1-combo mean: real = %.3f, rand = %.3f, ttest p = %.4f, ranksum p = %.4f\n', ...
    overall.real_1_mean, overall.rand_1_mean, overall.ttest_p_1, overall.ranksum_p_1);

fprintf('2-combo mean: real = %.3f, rand = %.3f, ttest p = %.4f, ranksum p = %.4f\n', ...
    overall.real_2_mean, overall.rand_2_mean, overall.ttest_p_2, overall.ranksum_p_2);

fprintf('3-combo mean: real = %.3f, rand = %.3f, ttest p = %.4f, ranksum p = %.4f\n', ...
    overall.real_3_mean, overall.rand_3_mean, overall.ttest_p_3, overall.ranksum_p_3);

%% Plot histograms
figure;
subplot(3,2,1); histogram(all_real_1, 'FaceAlpha', 0.6); title('Real z (nimg = 1)');
subplot(3,2,2); histogram(all_rand_1, 'FaceAlpha', 0.6); title('Rand z (nimg = 1)');
subplot(3,2,3); histogram(all_real_2, 'FaceAlpha', 0.6); title('Real z (nimg = 2)');
subplot(3,2,4); histogram(all_rand_2, 'FaceAlpha', 0.6); title('Rand z (nimg = 2)');
subplot(3,2,5); histogram(all_real_3, 'FaceAlpha', 0.6); title('Real z (nimg = 3)');
subplot(3,2,6); histogram(all_rand_3, 'FaceAlpha', 0.6); title('Rand z (nimg = 3)');

%% Summary table
T = struct2table(metrics);
disp(T)

% Plot bar comparisons
figure;
subplot(1,3,1); bar([T.real_1_mean, T.rand_1_mean]); title('nimg=1'); legend('Real','Rand'); xlabel('Drawing'); ylabel('z');
subplot(1,3,2); bar([T.real_2_mean, T.rand_2_mean]); title('nimg=2'); legend('Real','Rand'); xlabel('Drawing'); ylabel('z');
subplot(1,3,3); bar([T.real_3_mean, T.rand_3_mean]); title('nimg=3'); legend('Real','Rand'); xlabel('Drawing'); ylabel('z');


figure;
for i = 1:4
    subplot(2,2,i)
    histogram(all_top1_rand_per_drawing{i}, 'FaceColor', [0.4 0.6 1], 'EdgeColor', 'none');
    hold on;
    xline(all_top1_real_per_drawing(i), 'r', 'LineWidth', 2)
    title(sprintf('Drawing %d', i));
    xlabel('z-value'); ylabel('Count');
    legend('Rand top1/rep', 'top Real', 'Location', 'best')
end

figure;
drawings_to_plot = 5:8;
for i = 1:length(drawings_to_plot)
    d = drawings_to_plot(i); % actual drawing index
    subplot(2,2,i)
    histogram(all_top1_rand_per_drawing{d}, 'FaceColor', [0.4 0.6 1], 'EdgeColor', 'none');
    hold on;
    xline(all_top1_real_per_drawing(d), 'r', 'LineWidth', 2)
    title(sprintf('Drawing %d', d));
    xlabel('z-value'); ylabel('Count');
    legend('Rand top1/rep', 'Top Real', 'Location', 'best')
end
