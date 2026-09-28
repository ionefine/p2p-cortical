% Crop+Elche_Drawings.m
% Written by IF 3/7/2025
% Commened by ES 7/25/25
%
% takes in patient drawings and crops them


clear
addpath('../imregcode/');
%datadir = fullfile('Z:', 'IoneFine', 'MyProjects', 'p2p-cortical', 'ElcheDrawings', 'drawings');
datadir= fullfile(pwd);
dirList = dir([datadir, '/img*']); % finding the directories of patient drawings
nFiles = length(dirList); % number of patient drawings to crop



 for ff = 1:nFiles % for each folder of patient drawings

     % load original and edited images (you need both to run this program)
     orig_img = double(imread(fullfile(datadir, dirList(ff).name, 'drawings', [dirList(ff).name, '.png']))); % patient drawing
     edit_img = double(imread(fullfile(datadir, dirList(ff).name, 'drawings', [dirList(ff).name, '_edit.png']))); % patient drawing

     sz = size(orig_img); img = mean(orig_img, 3); % move to grayscale

     % find crop region for white space
     tmp = find(img(round(sz(1)/2), :)== 0);
     crop(2) = tmp(1)+1; crop(4) = tmp(end)-1;
     tmp = find(img(:, round(sz(1)/2)) == 0);
     crop(1) = tmp(1)+1; crop(3) = tmp(end)-1;

     orig_img_tmp = orig_img(crop(1):crop(3), crop(2):crop(4));
     edit_img_tmp = edit_img(crop(1):crop(3), crop(2):crop(4));

     % tighter crop
     [~, st_row, st_col, end_row, end_col] = crop_img(orig_img_tmp, 50);
     orig_img_crop = double(orig_img_tmp(st_row:end_row, st_col:end_col));
     edit_img_crop = double(edit_img_tmp(st_row:end_row, st_col:end_col));

     % save cropped version
     imwrite(orig_img_crop./max(orig_img_crop(:)), fullfile(datadir, dirList(ff).name,'drawings',  [dirList(ff).name, '_crop.png'])); % patient drawing
     imwrite(edit_img_crop./max(edit_img_crop(:)), fullfile(datadir, dirList(ff).name,'drawings',  [dirList(ff).name, '_edit_crop_0.png'])); % patient drawing

 end
