% Simulate_WFMA
%
% Simulation of indivual phosphenes using Wireless Floating Microelectrode
% Arrays (WFMAs)
%  Based on: Michael P Barry, Roksana Sadeghi, Vernon L Towle, Kelsey Stipp, Patricia Grant, Frank John Lane, Janet P Szlyk, Gislin Dagnelie, Philip R Troyk;
% Contributed Talk Session III: Characteristics of electrically-induced visual percepts in the first human with the Intracortical Visual Prosthesis. Journal of Vision, forthcoming 2023.
%
% Written by GMB & IF
% 25/02/2023 moved into clean folder (IF)
% August 2025 fixed plotting Eirini Schoinas

clear
clear all

rng(1171960)  % fix the random number generator. This affects the ocular dominance/orientation maps
c.I_k = 1000; % high electric field fall off -  only stimulating directly under the electrode

% define pulse train
tp = p2p_c.define_temporalparameters(); % define the temporal model
trl.amp = 60; trl.freq = 200;
trl.pw = 2*10^(-4);   trl.dur= .8;
trl = p2p_c.define_trial(tp,trl);

%  based on their abstract, the  sets of electrodes had the following visual field positions:
all_wfma(1).x = -3.5; all_wfma(1).y  = -3.5;
all_wfma(2).x = -4; all_wfma(2).y  = -8;
all_wfma(3).x = -30; all_wfma(3).y  = 0;

% sampling
v.pixperdeg =24;  %visual field map size and samping
c.pixpermm = 24; % default 6, resolution of electric field sampling, for very small electrodes may need to be decreased
scFac = 8;

% find out where the wfma are in cortex
for ii = 1
    if ii == 1
        c.cortexHeight = [0, 20]; % degrees top to bottom, degrees LR,
        c.cortexLength = [20, 55];
        v.visfieldHeight = [-10,0];
        v.visfieldWidth= [-10,0];
    elseif ii==2
        c.cortexHeight = [0, 20]; % degrees top to bottom, degrees LR,
        c.cortexLength = [20, 55];
        v.visfieldHeight = [-15,0];
        v.visfieldWidth= [-15,0];
    elseif ii==3
        c.cortexHeight = [-15, 15]; % degrees top to bottom, degrees LR,
        c.cortexLength = [20, 70];
        v.visfieldHeight = [-30,30];
        v.visfieldWidth= [-45,0];
    end

    v.e.x = all_wfma(ii).x;     v.e.y = all_wfma(ii).y;
    v = p2p_c.define_visualmap(v); % defines the visual map
    c = p2p_c.define_cortex(c); % define the properties of the cortical map
    [c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
    c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space
    
    % create and plot wfma
    wfma = create_WFMA(c.e.x,  c.e.y);
    plot_WFMA(wfma); % figure 100
    wfma = wfma([1, 2, end]); % pick first second and last electrode to simulate phosphenes
    
    % phosphenes for individual electrodes
    ct = 1;
    for e=1:3
        c.e.x = wfma(e).x; c.e.y = wfma(e).y;c.e.radius =  wfma(e).radius;
        v = p2p_c.define_visualmap(v); % defines the visual map
        c = p2p_c.define_cortex(c); % define the properties of the cortical map
        [c, v] = p2p_c.generate_corticalmap(c, v); % create ocular dominance/orientation/rf size maps on cortical surface
        v = p2p_c.c2v_define_electrodes(c,v); % convert electrode locations from cortex to visual space
        c = p2p_c.define_electrodes(c, v); % complete properties for each electrode in cortical space
        
        % set up the electrode locations in terms of their positions in the teeny array
        c = p2p_c.generate_ef(c); % generate map of the electric field for each electrode on cortical surface
                 figure(ii); subplot(1, length(wfma), e);
                 p2p_c.plotcortgrid(c.e.ef*2, c, gray(256), ii,['title(''electric field'')']); drawnow;
        %        savefig(['figures/Simulate_WFMA_Cortex_Fig', num2str(ii)]);

        % generate percepts for a pulse train
        v = p2p_c.generate_corticalelectricalresponse(c, v);  % create rf map for each electrode
        tmp_trl = p2p_c.generate_phosphene(v, tp, trl);
        figure(10+ii); subplot(1, length(wfma), ct)
        img = mean(tmp_trl.max_phosphene, 3);
        img = (img*scFac) +125;
        all_e_img(e, :,:)= img;

        % plotting individual phosphenes
        p2p_c.plotretgrid(img,  v, gray(256), 10+ii,['';]);
        a = gca; set(a, 'FontSize', 6);
        t = title(['E', num2str(e), 'LE' ]);    set(t, 'FontSize', 6);
        savefig(['figures/Simulate_WFMA_Visual Field_Fig', num2str(ii)]);
        ct = ct+1;
    end
end

% phosphenes for pairs  electrodes
figure(20 + ii); clf
subplot(1, 3, 1)
p2p_c.plotretgrid( squeeze(all_e_img(1, :, :)/2+all_e_img(2, :, :)/2), v, gray(256), 20+ii,['subplot(1, 3, 1)';]); hold on
set(t, 'FontSize', 6);    a = gca; set(a, 'FontSize', 6);
subplot(1, 3, 2)
p2p_c.plotretgrid( squeeze(all_e_img(1, :, :)/2+all_e_img(3, :, :)/2), v, gray(256), 20+ii,['subplot(1, 3, 2)';]); hold on
set(t, 'FontSize', 6);    a = gca; set(a, 'FontSize', 6);
subplot(1, 3, 3)
p2p_c.plotretgrid( squeeze(all_e_img(2, :, :)/2+all_e_img(3, :, :)/2), v, gray(256), 20+ii,['subplot(1, 3, 3)';]); hold on
set(t, 'FontSize', 6);    a = gca; set(a, 'FontSize', 6);
savefig(['figures/Simulate_WFMA_Visual FieldPair_Fig']);

%% WFMA functions
function plot_WFMA(wfma)
figure(100); clf
for i = 1:length(wfma)
    plot(wfma(i).x, wfma(i).y, 'ko', 'MarkerSize', wfma(i).radius*10^7); hold on
end
end

function e = create_WFMA(xoffset, yoffset)
%Troyk, P., Bredeson, S., Cogan, S., Romero-Ortega, M., Suh, S., Hu, Z., ... & Bak, M.
% % (2015, April). In-vivo tests of a 16-channel implantable wireless neural stimulator.
% % 2015 7th International IEEE/EMBS Conference on Neural Engineering (NER) (pp. 474-477). IEEE.
r1 = sqrt(500/pi)/1e+7; r2 = sqrt(1000/pi)/1e+7;
r3 = sqrt(1500/pi)/1e+7; r4 = sqrt(2000/pi)/1e+7;
rlist = [r1 r2 r3 NaN ...
    r3 r4 r1 r2 r4 ...
    r2 r3 r4 r1 ...
    r1 r2 r3 r4];

e2e = 0.4;

r_off = [.2 0 .2  .4];
ct = 1;
for r = [1:4] -1
    if r == 1
        cv = 0:4;
    else
        cv = 0:3;
    end
    for c = cv
        e(ct).y = (1.2-(r*0.4)) + yoffset;
        e(ct).x = (c*e2e)+r_off(r+1) + xoffset;
        e(ct).radius = rlist(ct);
        ct = ct+1;
    end
end
e = e([1:3, 5:17]); % remove the counter electrode
axis equal
end

