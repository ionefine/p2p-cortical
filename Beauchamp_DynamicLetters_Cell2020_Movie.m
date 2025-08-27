% Beauchamp_DynamicLetters_Cell 2020_Movie.m
%
% Replicates drawing task from Beauchamp dynamic stimulation paper
% Beauchamp, M. S., Oswalt, D., Sun, P., Foster, B. L., Magnotti, J. F., Niketeghad, S., ... & Yoshor, D. (2020).
% Dynamic stimulation of visual cortex produces form vision in sighted and blind humans. Cell, 181(4), 774-783.
% written IF & GMB
%
% loads in RF maps created by Beauchamp_DynamicLetters_Cell 2020_RF.m
% from a mat file Beauchamp_DynamicLetters_RF_Figure4.mat etc.
% 
% 25/02/2023 moved into clean folder (IF)
% 06/03/2023 Split RF generation and movie generation into separate files (IF) 
% code review, commenting + fixing save video/getframe so it works on any display by ES August 2025

% set up 
clear all
debugflag = 0;
rng(11171964)  % fix the random number generator. This affects the ocular dominance/orientation maps
fps = 30;

% experiment List
expList = { 'Figure 4-grid'}; %, 'Figure 6'}; % 'Figure 4'}; % grid is calculated location based on the array

% Load Beauchamp data
Te = readtable("datasets/Beauchamp_2020_data.xlsx", 'Sheet','ElectrodeLocations');
To = readtable("datasets/Beauchamp_2020_data.xlsx", 'Sheet','ElectrodeOrder');

%% Set up Visual and Cortical field

% visual field
v.visfieldHeight = [-25, 25]; 
v.visfieldWidth= [0,25]; 

if debugflag
    v.pixperdeg = 5; 
else
    v.pixperdeg = 12;
end
v = p2p_c.define_visualmap(v);

% cortex
c.cortexHeight = [-35,35]; % degrees top to bottom, degrees LR,
c.cortexLength = [-60, 5];

if debugflag
    c.pixpermm = 5; 
else
    c.pixpermm = 12; 
end
c.e.radius = 0.25;
c = p2p_c.define_cortex(c);

% generte map between them 
[c, v] = p2p_c.generate_corticalmap(c, v);

flag = 'Simultaneous'; % stimulation mode: 'Sequential' or 'Simultaneous'

for ex = 1:length(expList) % for each experiment 
    load(['datasets/Beauchamp_DynamicLetters_RF', expList{ex}]); % load RF maps

    % electrode locations and stimulation order for this experiment
    eid =strcmp(Te.experiment,expList{ex});
    Tloc = Te(eid, :);
    eid =strcmp(To.experiment,expList{ex});
    Torder = To(eid, :); % list of the letters/shapes simulated in this experiment

    % now for each letter, simulate phosphenes
    for l =1:size(Torder, 1) % for each letter

        % open the video file for that letter
        filename = ['movies/Beauchamp_Cell2020_',expList{ex}, Torder.letter{l}, flag];
        vid = VideoWriter([filename, '.avi']);
        vid.FrameRate = fps;    
        open(vid);

        % order of stimulation for that letter
        oList = str2num(Torder.order{l});
        img = zeros(size( saved(1).rfmap)); % zero out the frame
        trl.lag = .5; % initial lag
        
        % for each electrode,  as it is stimulated in order
        for o = 1:length(oList) 
            disp(['simulating electrode', num2str(o)]);
            disp(Tloc.electrode(oList(o)));
            
            % get trial parameters for this electrode
            trl.amp = Tloc.amp(oList(o));    
            trl.pw = Tloc.pw(oList(o));    
            trl.dur =  Tloc.dur(oList(o));    
            trl.freq = Tloc.freq(oList(o));

            % timing depends on mode
            if strcmp(flag, 'Sequential')
                trl.lag = 0.75 + ((o-1)*(Torder.pt_lag(l)+ Torder.pt_duration(l)));
            elseif strcmp(flag, 'Simultaneous')
                trl.lag = 0.75 ;
            end

            % total simulation duration
            trl.simdur  = max(trl.lag) + .75;

            % set up temporal model and trial
            tp = p2p_c.define_temporalparameters();
            trl = p2p_c.define_trial(tp,trl);

            tp.model = 'linear'; % weird hack to deal with the spatiotemporal nonlinearity
            trl = p2p_c.spike_model(tp, trl); % create response to the pulse        
            t_ds = round(linspace(1, length(trl.resp), trl.simdur*fps)); % downsample response to video rate

            % apply compression nonlinearity
            tp.model  = 'compression';
            for t = 1:length(t_ds)
                tmp =saved(oList(o)).rfmap*trl.resp(t_ds(t)); %  multiple space by time
                rfimg(o, t, :, :) = single(p2p_c.nonlinearity(tp, tmp)); % pass it through the response nonlinearity
            end
        end

        % now loop over video frames to save
        for t = 1:length(t_ds)

            % sum across electrodes to get final image at this frame
            img = squeeze(sum(double(rfimg(:, t, :, :)), 1));

            % plot visual space
            p2p_c.plotretgrid(img*500, v, gray(256), ex); 
            drawnow; 
            
            % get screensize so that each frame plots on your screen
            % this is needed for getframe() to work properly 
            scr = get(0,'ScreenSize');   % screen size in pixels
            figW = round(scr(3)*0.7);    % 70% of screen width
            figH = round(scr(4)*0.7);    % 70% of screen height
            figX = (scr(3)-figW)/2;      % center horizontally
            figY = (scr(4)-figH)/2;      % center vertically
            
            % set where frame displays
            set(gcf,'Units','pixels','Position',[figX figY figW figH]);

            % get frame and write to video
            frame = getframe(gca);
            writeVideo(vid,frame);
        end
        close(vid);
    end % finished that letter
end % finished that experiment

