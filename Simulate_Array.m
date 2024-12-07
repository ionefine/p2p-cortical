% Simulate_SmallElectrodes


% Written by GMB & IF
% 25/02/2023 clean up
% 12/7/2024 clean up

clear
clear all

% rng(1171960)  % fix the random number generator. This affects the ocular dominance/orientation maps

c.I_k = 1000;
% c.I_k = 1000 represents rapid fall off in the electric field as a function of distance from the electrode
%   i.e. only stimulating directly under the electrode.
% c.I_k = 6.75 is the default value. It represents fall off in the electric field as a function of distance based on Tehovnik 2006
% https://doi.org/10.1152/jn.00126.2006

% define pulse train
tp = p2p_c.define_temporalparameters(); % define the temporal model

trl.amp = 60; trl.freq = 200;
trl.pw = 2*10^(-4);   trl.dur= .8;
trl = p2p_c.define_trial(tp,trl);


% sampling of cortical and visual fields
v.pixperdeg = 24;  % visual field map samping
c.pixpermm = 24;   % resolution of electric field sampling

% Simulation of individual phosphenes using Wireless Floating Microelectrode Arrays (WFMAs)
%  Based on: Michael P Barry, Roksana Sadeghi, Vernon L Towle, Kelsey Stipp, Patricia Grant, Frank John Lane, Janet P Szlyk, Gislin Dagnelie, Philip R Troyk;
% Contributed Talk Session arrayI: Characteristics of electrically-induced visual percepts in the first human with the Intracortical Visual Prosthesis. Journal of Vision, forthcoming 2023.
%  Based on their abstract, the sets of electrodes had the following visual field positions:

% They had three possible arrays in different locations.
array = 1;

% define where the wfma is in cortex

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
    % v = p2p_c.define_visualmap(v); % defines the visual map
     c = p2p_c.define_cortex(c); % define the properties of the cortical map
    % 
    % [c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
     c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space

    array = create_array(c.e.x,  c.e.y);
    plot_array(array);
   

    for e = 1:length(array)
        c.e.x = array(e).x; c.e.y = array(e).y; c.e.radius =  array(e).radius;
        v = p2p_c.define_visualmap(v); % defines the visual map
        c = p2p_c.define_cortex(c); % define the properties of the cortical map
        [c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
        v = p2p_c.c2v_define_electrodes(c,v); % convert electrode locations from cortex to visual space
        c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space

        c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface
        
                figure(1); subplot(ceil(sqrt(length(array))), ceil(sqrt(length(array))), e);
                 p2p_c.plotcortgrid(c.e.ef*256, c, gray(256), array,['title(''electric field'')']); drawnow;
  
        % generate percepts 
        v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode
        tmp_trl = p2p_c.generate_phosphene(v, tp, trl);

        
        img = mean(tmp_trl.maxphos, 3); 
        img = img./max(abs(img(:))); % normalize so max is 1
        img = (img+.5)*127; %  scale 
        all_e_img(e, :,:) = img; % save the image in case you want to plot more than one electrode at a time
        
        figure(2); subplot(ceil(sqrt(length(array))), ceil(sqrt(length(array))), e)
        p2p_c.plotretgrid(img,  v, gray(256), 2,['';]);
        a = gca; set(a, 'FontSize', 6);
        t = title(['E', num2str(e), 'Phosphene' ]);    set(t, 'FontSize', 6);
    end


% calculat edistance between each pair of electrodes
for e1 = 1:length(array)
    for e2 =  1:length(array)
dd(e1, e2) = sqrt((array(e1).x-array(e2).x).^2 +(array(e1).y-array(e2).y).^2)
    end
end

% the predicted percept produced by pairs of electrodes
figure(3); clf
p2p_c.plotretgrid( squeeze(all_e_img(1, :, :)/2+all_e_img(2, :, :)/2), v, gray(256), 3,['subplot(1, 2, 1)';]); hold on
title(['e2e distance = ', num2str(dd(1,2)), ' mm']);
[i,j ] = find(dd == max(dd(:)), 1);
p2p_c.plotretgrid( squeeze(all_e_img(i, :, :)/2+all_e_img(j, :, :)/2), v, gray(256), 3,['subplot(1, 2,2)';]); 
title(['e2e distance = ', num2str(dd(1,end)), ' mm']);

function plot_array(wfma)
% plots the array configuration

figure(100); clf
for i = 1:length(wfma)
    plot(wfma(i).x, wfma(i).y, 'ko', 'MarkerSize', wfma(i).radius*10^7); hold on
end
end

function e = create_array(xoffset, yoffset)

% Here we are simulating the following array
% Troyk, P., Bredeson, S., Cogan, S., Romero-Ortega, M., Suh, S., Hu, Z., ... & Bak, M.
% % (2015, April). In-vivo tests of a 16-channel implantable wireless neural stimulator.
% % 2015 7th International IEEE/EMBS Conference on Neural Engineering (NER) (pp. 474-477). IEEE.
% 10.1109/NER.2015.7146662

% the array has varying sizes of electrode
r1 = sqrt(500/pi)/1e+7; r2 = sqrt(1000/pi)/1e+7;
r3 = sqrt(1500/pi)/1e+7; r4 = sqrt(2000/pi)/1e+7;
rlist = [r1 r2 r3 NaN ...
    r3 r4 r1 r2 r4 ...
    r2 r3 r4 r1 ...
    r1 r2 r3 r4];

e2e = 0.4; % electrode to electrode distances

r_off = [.2 0 .2  .4];
ct = 1;
for r = [0:3]
    if r == 1
        cv = 0:4;
    else
        cv = 0:3;
    end
    for c = cv
        e(ct).y = (1.2-(r*e2e)) + yoffset;
        e(ct).x = (c*e2e)+r_off(r+1) + xoffset;
        e(ct).radius = rlist(ct);
        ct = ct+1;
    end
end
e = e([1:3, 5:17]); % remove the counter electrode
end





