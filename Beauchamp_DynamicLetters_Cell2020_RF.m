% Beauchamp_DynamicLetters_Cell 2020.m
%
% Replicates drawing task from Beauchamp dynamic stimulation paper
% Beauchamp, M. S., Oswalt, D., Sun, P., Foster, B. L., Magnotti, J. F., Niketeghad, S., ... & Yoshor, D. (2020).
% Dynamic stimulation of visual cortex produces form vision in sighted and blind humans. Cell, 181(4), 774-783.
% written IF & GMB
%
%  creates RF maps and saves them mat file Beauchamp_DynamicLetters_RF_Figure4.mat etc.
%
% 25/02/2023 moved into clean folder (IF)
% 06/03/2023 Split RF generation and movie generation into separate files
% (IF)  and used electrode array instead of patient report to determine the
% location of phosphenes
% code review + commenting by ES August 2025

% set up 
clear
rng(11171964);  % fix the random number generator. This affects the ocular dominance/orientation maps

% Figure 4 and Figure 6 replicate the corresponding figures in Beauchamp's
% paper. Figure 4 - grid uses the fact that we know the electrode array has
% 2mm electrode separation and simulated the expected location of
% phosphenes, allowing a certain amount of mashing of the cortical surface
expList ={'Figure 4-grid'};% {'Figure 4', 'Figure 6', 


%% load Beauchamp data
Te = readtable("datasets/Beauchamp_2020_data.xlsx", 'Sheet','ElectrodeLocations');
To = readtable("datasets/Beauchamp_2020_data.xlsx", 'Sheet','ElectrodeOrder');

%% define cortical and visual space
% define visual field and sampling rate
v.visfieldHeight = [-35, 35]; 
v.visfieldWidth= [0,35]; 
v.pixperdeg = 18;
v = p2p_c.define_visualmap(v);

% define cortex
c.cortexHeight = [-35,35]; % degrees top to bottom, degrees LR,
c.cortexLength = [-80, 5];
c.pixpermm = 18; 
c.squish = 0.63;
c.a = .1474;
c = p2p_c.define_cortex(c);

% Generate baseline cortical feature maps
[c, v] = p2p_c.generate_corticalmap(c, v);

PLOT = 0;

%% find (and plot) the locations of the electrodes + find and save RF maps 
for ex = 1:length(expList) % for each experiment
    disp(['simulating ', expList{ex}]);
    subplot(1, length(expList), ex);

    % filter table rows for this experiment
    eid =strcmp(Te.experiment,expList{ex});
    Tloc = Te(eid, :);

    for e = 1:size(Tloc, 1) % for each electrode
        disp(['simulating electrode ', num2str(e)]);

        % write electrode parameters into v and c
        v.e.x = Tloc.x(e); v.e.y = Tloc.y(e);
        c.e.radius = Tloc.radius(e);

        % find electrode and electric field on cortex
        c = p2p_c.define_electrodes(c, v);
        c = p2p_c.generate_ef(c);

        % create RF map for each electrode
        v = p2p_c.generate_corticalelectricalresponse(c, v); 

        % save the rf map for that electrode, averaged across both eyes
        saved(e).rfmap = mean(v.e.rfmap, 3); 
        
        % if plotting 
        if PLOT
            p2p_c.plotretgrid(v.e.rfmap(:, :, 1)*255, v, gray(256), 1); drawnow
        end     

        % clean up to avoid errors
        rmfield(v,'e');
    end
    %% saving 
    % save RF maps for creating figures in
    % Beauchamp_DynamicLetters_Cell202_Movie
    save(['datasets/Beauchamp_DynamicLetters_RF', expList{ex}], 'saved');
    disp("saved")
end