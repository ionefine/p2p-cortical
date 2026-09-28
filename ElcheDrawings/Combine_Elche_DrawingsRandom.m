% Combine_Elche_Drawings.m
% Written by IF 3/7/2025
%
% takes in a patient drawing and finds a percept that looks similar

clear
close all

debugflag = 0;
addpath(genpath('../'));

if strcmp(computer, 'GLNXA64')
    datadir = fullfile(pwd);
else
    datadir = fullfile('Z:', 'IoneFine', 'MyProjects', 'p2p-cortical', 'ElcheDrawings');
end
cd(datadir)
dirList = dir([datadir, '/*img*']);
load(fullfile(datadir, 'patient_draw.mat'));
n_Drawings = length(patient_draw); % number of patient drawings to simulate
n_SimDrawCombined = 3; % maximum number of sites stimulated by a single electrode

%% do the fitting and matching


for d = 1:n_Drawings % for each patient folder drawing
    clear good_model
    disp(['working on drawing ', num2str(d), ' out of ', num2str(n_Drawings)]);
    clear tmp_model combos X B nimg corrval ssimval zval cmbx comboimg
    load(fullfile(datadir, dirList(d).name, 'models_8_1_2025', [dirList(d).name, '.mat']), 'sim_model');
   
    for rep = 1:50
       disp(['working on rep ', num2str(rep)]);
        load(fullfile(datadir, dirList(d).name, 'random_models', ['random_', dirList(d).name, '_', num2str(rep), '.mat']));

        % we don't care which subimage things came from, so unwrap and put them
        % in a single vector/matrix
        ct = 1; 
        for sd = 1:patient_draw(d).n_SubImages
            shuf_ind = randperm(ceil(patient_draw(d).nToSave/3));
            for ii = 1:ceil(patient_draw(d).nToSave/3)
                if rep==1
                good_model.img( :, ct) = double(sim_model.subimg{sd}(:,ii))/255;
                good_model.corr(ct) = sim_model.corr{sd}(ii);
                good_model.x(ct) = sim_model.x{sd}(ii);
                good_model.y(ct) = sim_model.y{sd}(ii);
                good_model.radius(ct) = sim_model.radius{sd}(ii);
                good_model.sub_draw_ind(ct) = sd-1;
                end
                good_model.randimg(:, ct) =double(sim_model_random.subimg{sd}(:,shuf_ind(ii)))/255;
                ct = ct+1;
            end
        end

        %% create all combos of up to n_SimDrawCombined images
        combos = dec2bin(1:(2^size(good_model.randimg, 2))-1)=='1'; % all combinations of model images
        ind = sum(combos, 2)<=n_SimDrawCombined;
        combos = combos(ind, :);
        clear AIC corrval zval nimg ssim r_cmbx cmbx comboimg r_nimg r_corrval r_ssimval r_zval
       
        parfor c = 1:size(combos, 1)
            if mod(c, 250)==0
                disp(['combination ', num2str(c), ' out of ', num2str(size(combos, 1))]);
            end
               X = good_model.img(:, combos(c,:));
            if rep==1
            B = X\patient_draw(d).subimg(1).img(:);
            nimg(c) = sum(combos(c,:)); % number of images in the combined model
            corrval(c) =corr(X*B, patient_draw(d).subimg(1).img(:));
            ssimval(c) = ssim(X*B, patient_draw(d).subimg(1).img(:), 'Exponents', [0 0 1]);
            zval(c) = atanh(corrval(c));
            cmbx(c, :) = combos(c,:);
            comboimg(c,:) = X*B;
            end
            Xr = good_model.randimg(:, combos(c,:));
            X(:, 1) =  Xr(:, 1);
            B = X\patient_draw(d).subimg(1).img(:);
            r_nimg(c) = sum(combos(c,:)); % number of images in the combined model
            r_corrval(c) =corr(X*B, patient_draw(d).subimg(1).img(:));
            r_ssimval(c) = ssim(X*B, patient_draw(d).subimg(1).img(:), 'Exponents', [0 0 1]);
            r_zval(c) = atanh(r_corrval(c));
            r_cmbx(c, :) = combos(c,:);
        end
        if rep==1
        [foo,id] = sort(zval,'descend');
        combo.nimg = nimg(id);
        combo.corrval = corrval(id);
        combo.zval = zval(id);
        combo.cmbx = cmbx(id,:);
        combo.comboimg = comboimg(id,:);
        end
        % re-sort based on correlation
        [foo,id] = sort(r_zval,'descend');
        r_combo(rep).nimg = r_nimg(id);
        r_combo(rep).corrval = r_corrval(id);
        r_combo(rep).zval = r_zval(id);
        r_combo(rep).cmbx = r_cmbx(id,:);
    end
    save(fullfile(datadir, dirList(d).name, 'best_combinations_8_4_2025.mat'), 'combo');
    save(fullfile(datadir, dirList(d).name, 'rand_combinations_8_4_2025.mat'), 'r_combo');
end
