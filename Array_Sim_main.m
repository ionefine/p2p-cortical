
function array = Array_Sim_main(array)
% Calls the functions needed to simulate various arrays and examine what
% the perceptual experience would be.
% send in a, must contain:
% array.spaceFac = 4 (or 3, 2, 1)
% array.arrayList =  {'optimal', ' regular_cortex','regular_visualfield'};%

% this function calls functions to generate electrode locations and RF maps
% based on the the array parameters
% to simulate what someone actually sees from a video run 
% Array_Sim_movie after

% written by IF and GMB
% code review, commented, cleaned up by Eirini Schoinas

%% set up

% change this to your home and data directories
params.homedir = '/Users/eschoinas/Desktop/p2p-cortical';
params.datadir = '/Users/eschoinas/Desktop/p2p-cortical/array';
params.plot = 1; % do you want plottting 
params.RF = 1; % do we want to generate the RF maps? or are they saved?
params.loc = 1; % do we need to generate electrode locations? or are they saved?


% up down then left right, we only generte a fourth of RF maps 
% and then mirror them to same time
params.flag_flip = [1, 1]; 

params.scale = [1 127]; % converts to uint8
params.overwrite = 1;

% set how many parral process you want to run with a parfor
% used when finding RF maps
% if par_size = 1 you dont use a parfor
par_size= 1;
disp(['using a parallel pool of ', num2str(par_size)]);

% different electrode set ups
if array.spaceFac == 4
    array.nelect = 79;
    c.pixpermm = 16; %16; % 16; % creates a rfmap of 936
    v.pixperdeg =16; %16; % 18;
elseif  array.spaceFac ==3
    array.nelect = 163;
    c.pixpermm = 16; %16; % 16; % creates a rfmap of 936
    v.pixperdeg =16; %16; % 18;
elseif array.spaceFac ==2
    array.nelect = 399;
    c.pixpermm = 14; %16; % 16; % creates a rfmap of 936
    v.pixperdeg =14; %16; % 18;
elseif array.spaceFac == 1
    array.nelect = 1880;
    c.pixpermm = 10; %16; % 16; % creates a rfmap of 936
    v.pixperdeg = 12; %16; % 18;
end

% eccentricty and cortex dimensions
array.eccLim = [0  32];
c.cortexHeight = [-40,40]; % degrees top to bottom, degrees LR,
c.cortexLength = [-80, 80];

% visual field dimensions
v.visfieldHeight = [-array.eccLim(2)*1.5,array.eccLim(2)*1.5];
v.visfieldWidth=[-array.eccLim(2)*1.5,array.eccLim(2)*1.5];

% size and elect electrical field fall off constant 
% c.I_k = 1000 represents rapid fall off in the electric field as a function of distance from the electrode
%   i.e. only stimulating directly under the electrode.
% c.I_k = 6.75 is the default value. It represents fall off in the electric field as a function of distance based on Tehovnik 2006
% https://doi.org/10.1152/jn.00126.2006
% use array.esize = small to deubug/run faster
if strcmp(array.esize, 'large')
    disp('simulating large electrodes')
    c.I_k =6.75;
    c.radius = 0.25;
else
    c.I_k = 1000;
    c.radius = 0.001;
end

% define cortex
c = p2p_c.define_cortex(c);

% start up parfor 
if par_size>1
    delete(gcp)
    parpool(par_size);
end

arrayList = array.arrayList; % array.arrayList = {'optimal','regular_cortex','regular_visualfield'};

% calcuate locations and RF maps for electrodes at different spacings

for aa = 1:length(array.arrayList)
    array.arrayStyle = arrayList{aa}; % hack so can't get overwritten
    array.filename = [params.datadir, filesep, array.arrayStyle, filesep, 'Array_Sim_', array.arrayStyle, '_', num2str(array.spaceFac), '_', array.esize];  
    
    %% calculate locations
    if params.loc == 1
        array = p2p_c.Array_Sim_Location(c, array);
        save([array.filename,'.mat'], 'array','-mat');
    else % load a previously generated file
        array = load([array.filename,'.mat'], '-.mat');
    end

    %% plot locations
    % we plot the electrodes the the visual field and cortex
    if params.plot == 1
        p.plotNums = [1, 3, 4, 2];
        p2p_c.Array_Sim_Plot(array, c, p);
        figure(3)
        set(gca, 'XLim', [-48 48]);
        set(gca, 'YLim', [-48 48]);
        figure(4)
        set(gca, 'XLim', [-5 80]);
        set(gca, 'YLim', [-30 30]);
    end

    %% generate RFs
   if params.RF == 1

       % if we want to compute for each electrode
        array.T = array.Tfull;

        % if we are flipping the array (because of symmtry of electrode locations)
        % then we may not need to compute for all the
        % electrodes
        if params.flag_flip(1)==1;  array.T  = array.T (array.T.vy>0, :); end %up down
        if params.flag_flip(2)==1;  array.T  = array.T (array.T.vx<0, :); % left right
        end
         
        % seperate electrodes for each parfor
        tmp = par_size*ceil(height(array.T)/par_size);
        range = reshape(1:tmp, tmp/par_size, par_size);

        % generate RF maps
        if par_size == 1 % one worker, no loop
            r = 1;
            rr = range(:, r);
            rfmaps = p2p_c.Array_Sim_GenerateRFmaps(array, v, c,params, rr);
            subset(r).rfmaps = rfmaps;
        else
            parfor r = 1:par_size % in parralel
                rr = range(:, r);
                disp(r);
                rfmaps = p2p_c.Array_Sim_GenerateRFmaps(array, v, c,params, rr);
                subset(r).rfmaps = rfmaps;
            end
        end

        % put subsets in parfors back together
        % put the maps back into the matrix
        array.rfmaps = subset(1).rfmaps;

        for r = 2:length(subset)
            array.rfmaps = cat(3, array.rfmaps, subset(r).rfmaps);
        end
        disp(['Finished RF creation ', array.arrayStyle]);
        save([array.filename, '_rfs.mat'], 'array', 'params', '-v7.3') % save to .mat file
   end
end

