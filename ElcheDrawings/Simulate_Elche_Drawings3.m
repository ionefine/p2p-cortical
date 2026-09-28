% Simulate_Elche_Drawings.msaved_models
% Written by IF 3/7/2025
% Commented by ES 7/25/25
%
% takes in a patient drawing and finds a percept that looks similar


clear
close all
addpath('../imregcode/');
addpath(genpath('../'));
if strcmp(computer, 'GLNXA64')
    datadir= fullfile(pwd);
else
    datadir = fullfile('Z:', 'IoneFine', 'MyProjects', 'p2p-cortical', 'ElcheDrawings');
end
cd(datadir)
dirList = dir([datadir, '/*img*']);
nDrawings = length(dirList); % number of patient drawings to simulate

n_SimDraw = 21;  % how many simulated drawings to save in total, divided equally among subimages

debugflag = 0;

flag_rectify = 1;
% two methods of fitting patient images
% one is to assume they only draw the bright bit of the images, when doing
% the fitting - though the final image shows bright and dark
% flag_rectify = 1, seems to work better
% the other doesn't rectify (flag_rectify = 0)

%% load patient drawings
sim_draw = [];
if nDrawings == 0
    error("failed to find patient drawing files");
end

for d = 1:nDrawings % for each patient folder drawing
    % sometimes we break images into pieces, how many subimages are there?
    subImgList = dir(fullfile(datadir, dirList(d).name,'drawings', '*edit_crop_*.png'));


    % store n_SubImages and nToSave for each drawing folder
    patient_draw(d).n_SubImages= length(subImgList);
    patient_draw(d).nToSave = round(n_SimDraw/patient_draw(d).n_SubImages);

    for sd = 1:patient_draw(d).n_SubImages % all the subdrawings
        % load subdrawing, create spatial refence object
        patient_draw(d).subimg(sd).img = double(imread(fullfile(datadir, dirList(d).name, 'drawings', subImgList(sd).name))); % patient drawing
        patient_draw(d).subimg(sd).fixed = imref2d(size(patient_draw(d).subimg(sd).img)); %reference for warping to patient_img

        % save space for model images
        sim_draw(d).subimg{sd} = uint8(rand([prod(size(patient_draw(d).subimg(sd).img)),n_SimDraw])*255);
        sim_draw(d).size{sd} = size(patient_draw(d).subimg(sd).img);
        sim_draw(d).corr{sd} = -ones(n_SimDraw,1);
        sim_draw(d).ssim{sd} = -ones(n_SimDraw,1);
        sim_draw(d).radius{sd} = NaN(n_SimDraw,1);
        sim_draw(d).x{sd} = NaN(n_SimDraw,1);
        sim_draw(d).y{sd} = NaN(n_SimDraw,1);
        sim_draw(d).subID{sd} = NaN(n_SimDraw,1);
    end
end

%% define the cortical model
tp = p2p_c.define_temporalparameters(); % define the temporal model

% because we are only interested in the shape of the percept, which doesn't
% depend on electrical stimulation parameters, we don't care what the
% exact pulse train is, so we set trl.freq to NaN which basically fixes
% response strength to 1, regardless of other stimulation params
trl.freq = NaN;
% trl.pw = 2*10^(-4);   trl.dur= .4;trl.amp = 60;
trl = p2p_c.define_trial(tp,trl);

% sampling of cortical and visual fields
v.pixperdeg = 15;  % visual field map samping, lowering this speeds things up
c.pixpermm = 15;   % resolution of electric field sampling

% cortical and visual field parameters
c.cortexHeight = [-20, 20]; % degrees top to bottom, degrees LR,
c.cortexLength = [20, 55];
c.onoff_ratio=.6;
v.visfieldHeight = [-6, 0];
v.visfieldWidth= [-6, 0];

v = p2p_c.define_visualmap(v); % defines the visual map

c = p2p_c.define_cortex(c); % define the properties of the cortical map
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface

% define the electrodes
ecc = 3.5; % select a random eccentricity, since we're not trying to predict size based on eccentricity
if debugflag
    nLoc = 3;
    r = 0.001;
else
    nLoc = 17; % number of x and y electrode locations
    r = [.005 .0075 .01 .25 .05 .1]; % different possible electrode radii
end
c.I_k = 1000;% fall off in the electric field
% define + store possible electrode locations + radii
theta = linspace(pi+pi/8, 3*pi/2-pi/8, nLoc);
[x, y] = pol2cart(theta, ecc);
[Ex, Ey, Er] = meshgrid(x, y, r);Ex = Ex(:); Ey = Ey(:); Er = Er(:);
% randomize stored order
ind = randperm(length(Ex));
Ex=Ex(ind); Ey = Ey(ind);Er = Er(ind);


for e = 1:length(Ex(:))
    disp(['working on electrode ', num2str(e), ' out of ', num2str(length(Ex(:)))]);
    c = p2p_c.define_cortex(c); % define the properties of the cortical map
    [c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface

    c.e.radius = Er(e);
    v.e.x = Ex(e); v.e.y = Ey(e);

    c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space
    c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface

    v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode
    trl = p2p_c.generate_phosphene(v, tp, trl); % generate phosphene for each electrode

    img = mean(trl.max_phosphene, 3); % model image
    img = crop_img(img, 20);
    img = img./max(abs(img(:))); % normalize so max is 1
    if flag_rectify
        img(img<0) = 0; % ignore dark bits
        img = img*255;
    else
        img = (img+.5)*127; %  scale
    end
    % at the end of this we have a phosphene produced an electrode, which
    % we call a 'model', next we compare these model to all the drawings.

    % what we need here is to 'best' 21 models
    % 'best' is defined in 2 ways: 1. having a strong correlation with the
    % drawing and 2. not being strongly correlated with other saved models.

    corrThresh = .95;
    clear  img_aligned
    parfor d = 1:nDrawings % for each patient folder drawing
        for sd = 1:length(patient_draw(d).subimg)
            [s, r] = findScaleRotationNGC_if(img,patient_draw(d).subimg(sd).img);
            [tform, peakcorr] = resolveSimilarityRotationAmbiguityNGC_if(img,patient_draw(d).subimg(sd).img, s, r);

            % is this model have a better correlation than the lowest so
            % far?  If so then either (1) replace the lowest in the list
            % with this one, or (2) find the correlation between this model
            % and those saved so far.  If the peak is > .95, then replace
            % this model with the highest correlated saved model.  Whee.

            if peakcorr > min(sim_draw(d).corr{sd})  % we have a model to save
                img_aligned = imwarp(img,tform,"OutputView",patient_draw(d).subimg(sd).fixed); % creates a model image in the patient image reference frame
                scaled_img = uint8(128*(img_aligned(:)-min(img_aligned(:)))/(max(img_aligned(:))-min(img_aligned(:))));
                if size(scaled_img) ~=size(patient_draw(d).subimg(sd).img(:))
                    error('model and patient images are different sizes');
                end
                % calculate the correlation between this model and the saved models
                allcorr = corr(double(sim_draw(d).subimg{sd}),img_aligned(:));
                if max(allcorr) > corrThresh % if too similar replace the highly correlated model
                    [foo,id] = max(allcorr);
                else % replace previous worst model with new one
                    [foo,id] = min(sim_draw(d).corr{sd});
                end

                % store models for now
                sim_draw(d).subimg{sd}(:,id) = scaled_img;
                sim_draw(d).radius{sd}(id) = c.e.radius;
                sim_draw(d).x{sd}(id) = v.e.x;
                sim_draw(d).y{sd}(id) = v.e.y;
                sim_draw(d).corr{sd}(id) = peakcorr;
                sim_draw(d).ssim{sd}(id) = ssim(img_aligned , patient_draw(d).subimg(sd).img, 'Exponents', [0 0 1]);
                sim_draw(d).subID = sd;
                sim_draw(d).size = size(img_aligned)
            end
        end
    end
end

if ~debugflag
    save(fullfile(datadir, 'patient_draw.mat'), 'patient_draw');
end

%% save models
for d = 1:nDrawings
    clear sim_model
    for sd = 1:patient_draw(d).n_SubImages% for each subimage
        % sort from best to worst
        [~,id] = sort(sim_draw(d).corr{sd},'descend');

        for i = 1:patient_draw(d).nToSave
            ii = id(i);
            sim_model.subimg{sd}(:, i) = [sim_draw(d).subimg{sd}(:,ii)];
            sim_model.radius{sd}(i) = [sim_draw(d).radius{sd}(ii)];
            sim_model.x{sd}(i) = [sim_draw(d).x{sd}(ii)];
            sim_model.y{sd}(i) = [sim_draw(d).y{sd}(ii)];
            sim_model.corr{sd}(i) = [sim_draw(d).corr{sd}(ii)];
            sim_model.ssim{sd}(i) = [sim_draw(d).ssim{sd}(ii)];
            sim_model.subID = sd;
            sim_model.size = sim_draw(d).size;
        end

    end
    % save models + model info in file
    save(fullfile(datadir, dirList(d).name, 'models', [dirList(d).name, '8_1_2025.mat']), 'sim_model');
end

