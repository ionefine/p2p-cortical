% Bosking_SizeAmplitude_JNeuro2017
%
% Simulates phosphene size as a function of electrode amplitude
%
% Saturation in Phosphene Size with Increasing Current Levels Delivered to Human Visual Cortex
% William H. Bosking, Ping Sun, Muge Ozker, Xiaomei Pei, Brett L. Foster, Michael S. Beauchamp, Daniel Yoshor
% Journal of Neuroscience 26 July 2017, 37 (30) 7188-7197; DOI: 10.1523/JNEUROSCI.2896-16.2017
%
% Written by GMB & IF
% 25/02/2023 moved into clean folder (IF)
% 03/03/2023 moved electrodes into table format and split into separate size vs.
% amplitude and size vs. eccentricity scripts (IF)
% 4/08/2023 commented ES

close all
n = 5;
Tloc = readtable('datasets/Bosking2017_data.xlsx'); % note some foveal electrodes are faked since couldn't be identified on the plot
eid = randperm(size(Tloc,1)); eid = eid(1:n); % randomly select n electrodes
Tloc = Tloc(eid,:);

%% define cortical and visual space
c.cortexHeight = [-15,15]; % degrees top to bottom, degrees LR,
c.cortexLength = [-80, 0];
c.pixpermm = 12;
c = p2p_c.define_cortex(c);

v.visfieldHeight = [-20,20]; v.visfieldWidth= [0,60]; v.pixperdeg = 12;
v = p2p_c.define_visualmap(v);
[c, v] = p2p_c.generate_corticalmap(c, v);

% temporal parameters
tp = p2p_c.define_temporalparameters();
tp.model = 'compression';
v.drawthr = 1;

% define parameters
amp = [50, 75, 100,150, 200, 250, 300, 350, 400];
dur =  200*10^-3*ones(size(amp));
pw =  .1 * 10^-3*ones(size(amp));
freq  = 200 *ones(size(amp));

% create table to hold parameter combos
Tsim = table(amp', pw', freq', dur');
Tsim.Properties.VariableNames = {'amp', 'pw','freq','dur'};

% set eccentricties (2 by default) + scaling factors 
e_loc = linspace(2, 12, 2); 
scList = exp(linspace(log(.2), log(.7), 7));

%% compute phosphene sized
for ii=1:length(e_loc) % loop over different eccentricities
    fprintf('Location %d of %d\n',ii,length(e_loc));

    for ss = 1:length(scList) % loop over different spatial spread values
        fprintf('Scale factor %d of %d\n',ss,length(scList));

        % Set electrode position in visual field
        v.e.ecc = e_loc(ii);
        v.e.ang = 0;

        % set spatial spread
        tp.sc_in = scList(ss);

        % generate trials with different amplitudes but fixed duration, freq, pw
        all_trl = p2p_c.loop_model(tp, Tsim);

        % Set electrode radius in mm
        c.e.radius = 0.25;

        % find electrode placement and electric field on the cortical map
        c = p2p_c.define_electrodes(c, v);
        c = p2p_c.generate_ef(c);

        % generate response
        v = p2p_c.generate_corticalelectricalresponse(c, v);

        for tt=1:length(amp) % loop over amplitudes

            trl = all_trl(tt); % grab trial    

            trl = p2p_c.generate_phosphene(v, tp, trl);
            img = mean(trl.max_phosphene, 3); % average over both eye-dominant cells

            % find phosphene size, brightness, maxresponse
            trl.sim_radius= mean([trl.ellipse(1).sigma_x trl.ellipse(1).sigma_y]);
            trl.sim_diameter = 2 * trl.sim_radius;
            trl.sim_brightness = max(trl.max_phosphene(:));
            trl.maxresp = max(trl.resp);   

            % store sizes for plotting          
            sim_sizes(ii,ss, tt) =  trl.sim_diameter;
            disp(['amp = ', num2str(amp(tt)), ' size ', num2str(sim_sizes(ii,ss,tt))]);

        end
    end
end

%% Plot normalized size as a function of Current (figure 3D)
disp('ready to plot')

% Line and symbol styles for different eccentricities, color, etc.
lineStyle = { '-', '--',};
symStyle = {'o', 's'};
colorList = parula(length(scList)-2);
figure(2); clf; hold on

% loop over eccentricity conditions
for ii = 1:length(e_loc)
    for ss = 1:length(scList)-4 % loop oversubset of scale factors
        max_size = max(sim_sizes(ii, ss,  :)); % take the max over all the amplitudes
        y = squeeze(sim_sizes(ii, ss, :)./max_size); % normalize so max is 1

        % plot points and line
        plot(amp,y ,'LineWidth',1, 'LineStyle', lineStyle{ii}, 'Color', colorList(ss,:));
        plot(amp,y,symStyle{ii},'Color',colorList(ss, :),'MarkerFaceColor',colorList(ss, :),'MarkerEdgeColor',colorList(ss, :),'MarkerSize',11,'LineWidth',1); hold on
    end
end

% labels, title, etc.
xlabel('Current (mA)');
ylabel('Normalized phosphene size');
set(gca,'XLim',[0,400]);
title('Normalized phosphene size vs current');

%set(gca,'YLim',[0 1]);
set(gca,'FontSize',8);


