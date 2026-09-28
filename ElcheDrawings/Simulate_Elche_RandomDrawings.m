
function Simulate_Elche_RandomDrawings(replist)
% Simulate_Elche_Drawings.msaved_models
% Written by IF 3/7/2025
% takes in a patient drawing and finds a percept that looks similar


addpath('../imregcode/');
addpath(genpath('../'));
if strcmp(computer, 'GLNXA64')
    datadir= fullfile(pwd);
else
    datadir = fullfile('Z:', 'IoneFine', 'MyProjects', 'p2p-cortical', 'ElcheDrawings');
end
cd(datadir)
dirList = dir([datadir, '/*img*']);
load(fullfile(datadir, 'patient_draw.mat'));
n_Drawings = length(patient_draw);
n_Reps = 100;

%% define cortex
tp = p2p_c.define_temporalparameters(); % define the temporal model
trl.freq = NaN;
trl = p2p_c.define_trial(tp,trl);
master_trl = trl;

% sampling of cortical and visual fields
v.pixperdeg = 15;  % visual field map samping, lowering this speeds things up
c.pixpermm = 15;   % resolution of electric field sampling

% define visual and cortical field dimensions
c.cortexHeight = [-20, 20]; % degrees top to bottom, degrees LR,
c.cortexLength = [20, 55];
c.onoff_ratio=.6;
v.visfieldHeight = [-6, 0];
v.visfieldWidth= [-6, 0];

v = p2p_c.define_visualmap(v); % defines the visual map
c = p2p_c.define_cortex(c); % define the properties of the cortical map
c.I_k = 1000;
v_master = v;
flag_rectify = 1;

for rep = 1:100
    disp(['working on rep ', num2str(rep), ' out of ', num2str(n_Reps)])
      [c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
       for d = 1:n_Drawings
        disp(['working on drawing ', num2str(d), ' out of ', num2str(n_Drawings)])
        load(fullfile(datadir, dirList(d).name, 'models', [dirList(d).name, '.mat']), 'sim_model');
        pd = patient_draw(d);
        sim_model_random = sim_model;
        for sd = 1:patient_draw(d).n_SubImages
            for i = 1:patient_draw(d).nToSave        
                c.e.radius = sim_model_random.radius{sd}(i);
                v.e.x = sim_model_random.x{sd}(i); v.e.y = sim_model_random.y{sd}(i);
                c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space
                c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface
                v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode
                trl = p2p_c.generate_phosphene(v, tp, trl);

                img = mean(trl.max_phosphene, 3); % model image
                img = crop_img(img, 20);
                img = img./max(abs(img(:))); % normalize so max is 1
                if flag_rectify
                    img(img<0) = 0; % ignore dark bits
                    img = img*255;
                else
                    img = (img+.5)*127; %  scale
                end

                [s, r] = findScaleRotationNGC_if(img,pd.subimg(sd).img);
                [tform, peakcorr] = resolveSimilarityRotationAmbiguityNGC_if(img,pd.subimg(sd).img, s, r);
                img_aligned = imwarp(img,tform,"OutputView",pd.subimg(sd).fixed); % creates a model image in the patient image reference frame
                scaled_img = uint8(128*(img_aligned(:)-min(img_aligned(:)))/(max(img_aligned(:))-min(img_aligned(:))));
                if size(scaled_img) ~=size(pd.subimg(sd).img(:))
                    error('model and patient images are different sizes');
                end

                sim_model_random.subimg{sd}(:,i) = scaled_img;
                sim_model_random.corr{sd}(i) = peakcorr;
                sim_model_random.ssim{sd}(i) = ssim(img_aligned , pd.subimg(sd).img, 'Exponents', [0 0 1]);
            end
        end
        save(fullfile(datadir, dirList(d).name, 'random_models', ['random', dirList(d).name, '_', num2str(rep), '.mat']), 'sim_model_random');
    end
end


