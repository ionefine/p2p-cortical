% Simulate_SmallElectrodes


% Written by GMB & IF
% 25/02/2023 clean up
% 12/7/2024 clean up
% 6/24/2025 clean up (ES)

% clear workspace so its clean 
clear
clear all



% fix the random number generator. This affects the ocular dominance/orientation maps
rng(1171960)

% set electrical field fall off constant 
% c.I_k = 1000 represents rapid fall off in the electric field as a function of distance from the electrode
%   i.e. only stimulating directly under the electrode.
% c.I_k = 6.75 is the default value. It represents fall off in the electric field as a function of distance based on Tehovnik 2006
% https://doi.org/10.1152/jn.00126.2006
% c.I_k = 1000;
% c.I_k = 6.75;
c.I_k = 1000;

% define pulse train
tp = p2p_c.define_temporalparameters(); % define the temporal model

% define trial parameters
trl.amp = 60; trl.freq = 200; % amplitude, frequency
trl.pw = 2*10^(-4); trl.dur= .8; %pulse width, pulse duration 

% define entire pulse train waveform for trial
trl = p2p_c.define_trial(tp,trl);  

% sampling of cortical and visual fields
v.pixperdeg = 24;  % visual field map sampling
c.pixpermm = 24;   % resolution of electric field sampling

% Simulation of individual phosphenes using Wireless Floating Microelectrode Arrays (WFMAs)
% Based on: Michael P Barry, Roksana Sadeghi, Vernon L Towle, Kelsey Stipp, Patricia Grant, Frank John Lane, Janet P Szlyk, Gislin Dagnelie, Philip R Troyk;
% Contributed Talk Session arrayI: Characteristics of electrically-induced visual percepts in the first human with the Intracortical Visual Prosthesis. Journal of Vision, forthcoming 2023.
% Based on their abstract, the sets of electrodes had the following visual field positions:

% They had three possible arrays in different locations.
% define where the wfma is in cortex
array = 1;

if array == 1
    c.cortexHeight = [0, 20]; % degrees top to bottom, degrees LR,
    c.cortexLength = [20, 55];
    v.visfieldHeight = [-10,0];
    v.visfieldWidth= [-10,0];
    v.e.x = -3.5; v.e.y = -3.5;
elseif array==2
    c.cortexHeight = [0, 20]; % degrees top to bottom, degrees LR,
    c.cortexLength = [20, 55];
    v.visfieldHeight = [-15,0];
    v.visfieldWidth= [-15,0];
    v.e.x = -4; v.e.y = -8;
elseif array==3
    c.cortexHeight = [-15, 15]; % degrees top to bottom, degrees LR,
    c.cortexLength = [20, 70];
    v.visfieldHeight = [-30,30];
    v.visfieldWidth= [-45,0];
    v.e.x = -30; v.e.y = 0;
end


% define visual and coritcal maps
c = p2p_c.define_cortex(c); % define the properties of the cortical map
c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space
v = p2p_c.define_visualmap(v); % defines the visual map
c = p2p_c.define_cortex(c); % define the properties of the cortical map
[c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface

   
% create electrode array 
array = create_array(c.e.x,  c.e.y);
plot_array(array); % plot array 

% allocate array to store generated phosphenes for each electrode image and
% electric field image
all_e_img = uint8(ones(size(v.X,1), size(v.X, 2), length(array)));
all_e_elec = uint8(ones(size(c.X,1), size(c.X, 2), length(array)));

% all_e_elec = cell(1, length(array));


c.e = array; % list of electrodes

v = p2p_c.c2v_define_electrodes(c,v); % convert electrode locations from cortex to visual space
c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space
c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface

v = p2p_c.generate_corticalelectricalresponse_vectorized(c, v);  % create rf map for each electrode


trl_array = p2p_c.generate_phosphene(v, tp, trl); % generate phosphene for electrode

for e = 1:length(array)
    
    img = mean(trl_array(e).max_phosphene, 3);

    img = img./max(abs(img(:))); % normalize so max is 1
    img = (img+.5)*127; %  scale
    all_e_img(:,:, e) = uint8(img); % save the image in case you want to plot more than one electrode at a time
    all_e_elec(:,:, e) = uint8(c.e(e).ef * 256);
end


% plot subeset of images in grid
plot_img_subset(all_e_elec, 16, 'First 16 Electric Fields');
plot_img_subset(all_e_img, 16, 'First 16 Phosphenes');

% save video file of image arrays
plot_movie(all_e_elec, 'electric.avi')
plot_movie(all_e_img, 'phosphenes.avi')
disp(['Temporal setup completed in ', num2str(toc(section_tic), '%.2f'), ' seconds.']);

save('RFs4ML');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions %%

function plot_array(wfma)
    % plots the array configuration 
    figure(200); clf
    for i = 1:length(wfma)
        plot(wfma(i).x, wfma(i).y, 'ko', 'MarkerSize', wfma(i).radius*1e5); hold on
    end
    axis equal
    axis tight
end

function e = create_array(xoffset, yoffset)
% Here we are simulating a high density array with small electrodes
% the array has varying sizes of electrode
% 21 x 21 = 441 electrodes
e2e = 0.5; % electrode to electrode distances
ct = 1;
for r = -10:10
    for c = -10:10
        e(ct).y = (r*e2e) + yoffset;
        e(ct).x = (c*e2e) + xoffset;
        e(ct).radius = .0001;
        ct = ct+1;
    end
end
end

function e = plot_img_subset(all_img, n_plot, ttl)
    % plot subest of 3D phosphene image array
    % all_img: array images to plot
    % n_plot: number of images to plot
    % ttl: overall plot title
    
    figure;
    for e = 1:n_plot
        subplot(4, 4, e);
        imshow(all_img(:,:,e), []);
        title(['Electrode ' num2str(e)]);
    end
    sgtitle(ttl); 
end

function plot_movie(all_img, filename, framerate)
    % generates a movie from img array
    % all_img: 3D array of imgs (X × Y × num_electrodes)
    % filename: name of the output video file
    % framerate: optional, frames per second

    if nargin < 3
        framerate = 5;
    end

    v = VideoWriter(filename);
    v.FrameRate = framerate;
    open(v);

    for e = 1:size(all_img, 3)
        frame = all_img(:, :, e);
        rgb = repmat(mat2gray(frame), [1 1 3]);  % ensure grayscale & RGB
        writeVideo(v, rgb);
    end

    close(v);
    disp(['Saved movie to ', filename]);
end







