% CompareRandomRealCorrelations.m
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
load([fullfile(datadir, 'patient_draw.mat')]);
n_Drawings = length(patient_draw);
%% do the fitting and matching

for d = 1:n_Drawings
    load([fullfile(datadir, dirList(d).name, "best_combinations.mat")]);
    load([fullfile(datadir, dirList(d).name, "rand_combinations.mat")]);
    topnum = length(combo(1).zval);
    zval2 = NaN(length(r_combo), topnum);
    zval3 = NaN(length(r_combo), topnum);

    for rep = 1:length(r_combo)

        ind =find(r_combo(rep).nimg(1:topnum)==2);
        if ~isempty(ind)
            zval2(rep,1:length(ind)) = r_combo(rep).zval(ind);
        end

        ind =find(r_combo(rep).nimg(1:topnum)==3);
        if ~isempty(ind)
            zval3(rep,1:length(ind)) = r_combo(rep).zval(ind);
        end

    end
    figure(d)
    subplot(1,2,1)
    histogram(combo(1).zval(1:topnum), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.2);

    %hist(combo(1).zval(1:topnum));
    subplot(1,2,2)
    histogram(zval2(:), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.2, 'DisplayName', '2-image');
    %hist(zval2(:), zval3(:));
    drawnow;



end