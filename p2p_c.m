% p2p_c
%
%  holds all support functions for p2p cortex
% project.
%
% functions can be called from outside with 'p2p_c.<function name>'
% written IF and GMB
%
% 25/02/2023 moved into clean folder (IF)


classdef p2p_c
    methods(Static)
        %% definitions - temporal properties

        function trl = define_trial(tp, varargin)
            % creates a pulse train for a given trial
            % takes as input:
            %  tp - the temporal parameters of the model (at a minimum
            %  needs tp.dt, the temporal sampling)
            %  [trl] if you want to preset some trial values
            %  If you don't want to use a temporal model, set the freq to
            %  NaN
            %
            % commented 12/7/2024 IF
            % commented 6/24/2024 ES

            % if no trl struct passed create one, else use passed struct
            if nargin < 2;  trl = [];
            else; trl = varargin{1}; end
            
            % set default electrode index
            if ~isfield(trl,'e'); trl.e = 1; end

            % set electrical stimulation and simulation durations
            % durations are all in seconds
            if ~isfield(trl, 'dur');    trl.dur = 1000*10^-3;   end % electrical stimulation duration
            if ~isfield(trl, 'simdur');    trl.simdur = 3 ;   end % simulation duration, needs to be
            % longer to allow time for the neural response

            % construct time vector
            trl.t = 0:tp.dt:trl.dur-tp.dt;

            % define pulse parameters
            if ~isfield(trl, 'pw');     trl.pw = .1 * 10^-3;    end % pulse width 
            if ~isfield(trl, 'ip');     trl.ip = 0;             end % interphase delay
            if ~isfield(trl, 'lag');    trl.lag = 2*trl.pw;     end % delay before the pulse train begins
            if ~isfield(trl, 'order');  trl.order = 1;          end % 1 = cathodic first, -1  = anodic first
            if ~isfield(trl, 'freq');   trl.freq = 60;          end % NaN if not using a temporal model
            if ~isfield(trl, 'amp');    trl.amp = 100;          end % current amplitude in microAmps

            % generate pulse train 
            trl = p2p_c.generate_pt(trl, tp);

            % calculate charges
            trl.CperTrial = (trl.amp/1000) * trl.dur * trl.freq * trl.pw*10.^3; % charge per trial
            trl.CperPulse = trl.pw * trl.amp/1000; % charge per pulse
        end
        
        function trl = generate_pt(trl, tp)
            % generate pulse train given trial parameters 
            % and temporal model parameters
            % commented 6/24/2024 ES
            
            % if 'trl.on' and 'trl.off' specified convert to lag and duration
            if isfield(trl, 'on')
                trl.lag = trl.on;
                trl.dur = trl.off-trl.on;
            end

            if isnan(trl.freq)
                % if all you are interested in is space then don't use a
                % temporal model at all, space and time are separable and
                % it's much faster
                trl.pt = 1;
            else
                if trl.ip == 0  % no interphase delay 
                    on =  mod(trl.t,1/trl.freq) <=(trl.pw*2); % turn it on on
                    off = mod(trl.t-trl.pw,1/trl.freq) <=trl.pw & on;  % '& on' hack added by gmb;
                    tmp  = trl.amp.*(on-(2*off));
                else 
                    on =  mod(trl.t,1/trl.freq) <trl.pw;
                    delay =  (trl.pw+trl.ip); % time difference between on and off
                    off = mod(trl.t-delay,1/trl.freq) < trl.pw;
                    tmp  = trl.amp.*(on-off);
                end

                lag = round(trl.lag/tp.dt); % delay before beginning the pulse train (frames)
                trl.pt= zeros(1, lag+length(tmp));
                trl.pt(lag+1:lag+length(tmp))=tmp;

                trl.t = 0:tp.dt:(trl.dur+trl.lag); % include the lag
                trl.t = trl.t(1:end-1);

            end
            if trl.dur<trl.simdur % usually we simulate a little longer than the trial, to allow for the response
                trl.pt((end+1):round((trl.simdur/tp.dt))) = 0;
                trl.t = 0:tp.dt:trl.simdur-tp.dt;
            end
        end

        function c = define_cortex(c)
            % defines the parameters of the cortical sheet and the models
            % that define what the receptive fields will look like
            %
            % commented 12/7/2024 IF

            % typical log z transformation parameters (based on early
            % Schwartz model
            if ~isfield(c, 'efthr'); c.efthr = 0.05; end % electric field values below this are assumed to be zero
            if ~isfield(c, 'animal') ;  c.animal = 'human'; end
            if strcmp(c.animal, 'human')
                c.k = 15; %imcale
                if ~isfield(c, 'a'); c.a = 0.5; end %fovea expansion for human, macaque is 0.3
                c.shift = c.k*log(c.a);
                if ~isfield(c, 'squish');   c.squish = 1;  end % some cortices are just a little rounder than others, no judgment
                if ~isfield(c, 'cortexHeight'); c.cortexHeight = [-40,40]; end %[height in mm of cortex, 0 is midline)
                if ~isfield(c, 'cortexLength'); c.cortexLength = [-5, 80]; end %[length in mm of cortex, 0 is fovea)
                if ~isfield(c, 'pixpermm');  c.pixpermm = 8; end    % choose the resolution to sample in mm.
            elseif strcmp(c.animal, 'macaque')
                c.k = 5; c.squish = 1; %scale
                c.a = 0.3; % values set by eyeballing Toottell data
                c.shift = c.k*log(c.a);
                if ~isfield(c, 'cortexHeight'); c.cortexHeight = [-20,20]; end %[height in mm of cortex, 0 is midline)
                if ~isfield(c, 'cortexLength'); c.cortexLength = [-5,30]; end %[length in mm of cortex, 0 is fovea)
                if ~isfield(c, 'pixpermm');  c.pixpermm = 8; end    % choose the resolution to sample in mm.
            elseif strcmp(c.animal, 'mouse')
                errordlg('Sorry no model for mouse yet');
                %c.k = 1/40; % scale Garrett, 2014 FOR V1 how many mm of cortex represents 1 degree of visual field
            end

            % size and structure of receptive fields

            if ~isfield(c, 'rfmodel');   c.rfmodel = 'ringach';  end
            if ~isfield(c, 'rfsizemodel')
                c.rfsizemodel = 'keliris';  %  estimate of how sizes of rfs change as a function of eccentricity Keliris et al. PNAS 2019, 10.1073/pnas.1809612116
            end
            if strcmp(c.animal, 'human')
                if strcmp(c.rfsizemodel, 'keliris') % using Keliris electrophysiology from supplementary table 1
                    c.slope = 0.08; % 0.05; % in terms of sigma of a Gaussian
                    c.intercept = 0.16; % 0.69;
                    c.min = 0;
                elseif strcmp(c.rfsizemodel, 'bosking') %  Saturation in Phosphene Size with Increasing Current Levels Delivered to Human Visual Cortex, Bosking et al. J Neurosci 2017
                    c.slope =   .2620/2; % Bosking data is in terms of diameter, so take the values from Figure 4 (slope = 0.2620 and intercept  = 0.1787) and divide by 2
                    c.intercept = .1787/2;
                    c.min = 0;
                elseif strcmp(c.rfsizemodel, 'winawer')
                    c.slope = .1667;
                    c.min = 1.11;
                    c.intercept =  0.0721;
                end
            elseif strcmp(c.animal, 'macaque')
                if strcmp(c.rfsizemodel, 'keliris')
                    c.slope =  0.08; % in terms of sigma
                    c.intercept = 0.16;
                    c.min = 0;
                else
                    c.slope =  0.06; % in terms of sigma
                    c.intercept = 0.42;
                    c.min = 0;
                end
            elseif strcmp(c.animal, 'mouse') %
                errordlg('mouse model isn''t defined')
                c.intercept = 20;  % Check this ezgi
            end

            %% size and structure of receptive fields
            if ~isfield(c,'ar'); c.ar = 0.25; end % aspect ratio elongated rfs, based on Ringach 2002, J. Neurophysiology, 10.1152/jn.2002.88.1.455

            if ~isfield(c, 'delta'); c.delta = 2; end % this describes the distribution describing the separation between on and off receptive fields, Based Mata & Ringach, 2005 10.1152/jn.00668.2004

            if ~isfield(c, 'onoff_ratio'); c.onoff_ratio  = 0.8; end % off cells contribute less to perception than on cells

            % Ocular dominance structure based on Adams et al. 2007, 10.1523/JNEUROSCI.2923-07.2007
            if ~isfield(c, 'sig'); c.sig = .5;  end % The larger sig, the more the distribution of ocular dominance columns tends toward being 0 or 1
            if strcmp(c.animal, 'human')
                c.ODsize = 0.863; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
                c.filtSz = 3; % 3mm creates the initial OD and orientation maps
            elseif strcmp(c.animal, 'macaque')
                c.ODsize = 0.531; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
                c.filtSz = 1.85; % 3mm creates the initial OD and orientation maps
            elseif strcmp(c.animal, 'mouse')
                c.ODsize = NaN; % Adams 2007 Complete pattern of ocular dominance columns in human primary visual cortex, average width of a column in mm
                c.filtSz = NaN; % 3mm creates the initial OD and orientation maps
            end

            c.gridColor = [1,1,0]; % for drawing lines on cortex.
        end

        function [c] = generate_ef(c, varargin)
            % generates an electric field for each electrode on the
            % cortical surface (can take multiple at a time as input)
            % max electric field is normalized to 1 at this point
            %
            % commented 12/7/2024 IF
            % commented 6/26/2025 ES

            idx = 1:length(c.e);
          
            % set default electric field model if not specified
            if ~isfield(c, 'emodel')
                c.emodel = 'Tehovnik';      end

            for ii=1:length(idx)
                % compute distance matrix R from electrode center to each cortical grid point
                R=sqrt((c.X-c.e(idx(ii)).x).^2+(c.Y-c.e(idx(ii)).y).^2);

                % subtract radius to get effective distance; clamp to 0 when within electrode
                Rd = R-c.e(idx(ii)).radius; Rd(Rd<0) = 0;

                % initialize electric field
                pt_ef = ones(size(c.X));

                % compute electric field based on selected model
                if strcmp(c.emodel, 'Tehovnik')
                    % https://doi.org/10.1152/jn.00126.2006
                    if ~isfield(c, 'I_0');     c.I_0  = 1; end
                    if ~isfield(c, 'I_k');     c.I_k  = 6.75; end % controls rate of decay of current spread in cm
                    pt_ef=c.I_0./(1+c.I_k*Rd.^2);
                else
                    errordlg('electric field model not specified in code');
                end
                c.e(idx(ii)).ef = uint8(255.*pt_ef./max(pt_ef(:))); % uint8 to save space when making movies
            end
        end

        % generate functions
        function [c, v] = generate_corticalmap(c, v)
            % creates cortical maps for:
            %       ocular dominance (c.ODmap)
            %       orientation (c.ORmap)
            %       rf size (c.RFsizemap)
            %       on vs. off (c.ONOFFmap)
            %       distance between on and off subunits, small values represent complex cells, larger represent simple (c.DISTmap)
            %
            % commented 12/7/2024 IF

            % define cortex meshgrid
            c.x = linspace(min(c.cortexLength),max(c.cortexLength), (max(c.cortexLength)-min(c.cortexLength))*c.pixpermm);
            c.y = linspace(min(c.cortexHeight),max(c.cortexHeight), (max(c.cortexHeight)-min(c.cortexHeight))*c.pixpermm);
            [c.X,c.Y] = meshgrid(c.x,c.y);
            sz = size(c.X);

            %% Make the orientation and OD maps by bandpassing random noise

            % Rojer and Schwartz' method of bandpassing random noise:
            % Rojer, A.S. and E.L. Schwartz, Cat and monkey cortical columnar patterns
            %modeled by bandpass-filtered 2D white noise. Biol Cybern, 1990. 62(5): c. 381-91.
            %Make random noise: complex numbers where the angle is the orientation
            Z = exp(sqrt(-1)*rand(sz)*pi*2);

            % filter the noise to create initial columns
            freq = 1/c.ODsize; %cycles/mm (Try zero for big columns)
            filtPix = ceil(c.filtSz*c.pixpermm);
            [X,Y] = meshgrid(linspace(-c.filtSz/2,c.filtSz/2,filtPix),linspace(-c.filtSz/2,c.filtSz/2,filtPix));
            R = sqrt(X.^2+Y.^2);
            FILT = exp(-R.^2/c.sig.^2).*cos(2*pi*freq*R);  %Gabor

            %Convolve z with the filter
            W = conv2(Z,FILT,'same');
            c.ORmap = angle(W);
            WX = gradient(W);
            Gx = angle(WX);
            c.ODmap = normcdf(Gx*c.sig);


            %% Make the on and off  maps by bandpassing random noise
            % The idea that these vary smoothly is based on
            % Najafian, S., Koch, E., Teh, K.L. et al. A theory of cortical map formation in the visual brain.
            % Nat Commun 13, 2303 (2022). https://doi.org/10.1038/s41467-022-29433-y
            % in this paper the maps are related to orientation and ocular
            % dominance, but because that doesn't matter for our
            % simulations we're just generating new maps

            % filter the noise to create initial columns
            freq = 1/c.ODsize *2; %cycles/mm, doubling the frequency
            filtPix = ceil(c.filtSz/2*c.pixpermm);
            [X,Y] = meshgrid(linspace(-c.filtSz/2,c.filtSz/2,filtPix),linspace(-c.filtSz/2,c.filtSz/2,filtPix));
            R = sqrt(X.^2+Y.^2);
            FILT = exp(-R.^2/c.sig.^2).*cos(2*pi*freq*R);  %Gabor

            %Convolve z with the filter
            W = conv2(Z,FILT,'same');
            u = (angle(W)/(pi)); % distance map
            tmp = zeros(size(u));
            tmp(u>0) = -log(u(u>0))/c.delta; % exponential fall of of d, as described by Mata & Ringach 2005
            tmp(u<0) = log(-u(u<0))/c.delta;
            c.DISTmap = tmp;
            WX = gradient(W);
            Gx = angle(WX);
            c.ONOFFmap = normcdf(Gx*c.sig);

            % create angle and eccentricity maps
            [c.v.X ,  c.v.Y] = p2p_c.c2v_real(c, c.X, c.Y);
            [c.v.ANG, c.v.ECC] = cart2pol(c.v.X,c.v.Y);
            c.v.ANG = c.v.ANG *180/pi;

            % create a mesh
            v.zAng = linspace(0,max(v.eccList),v.n)'*exp(sqrt(-1)*v.angList*pi/180);
            c.v.gridAng = p2p_c.v2c_cplx(c, v.zAng);

            v.zEcc = (v.eccList'*exp(sqrt(-1)*linspace(-90,90,v.n)*pi/180))';
            c.v.gridEcc = p2p_c.v2c_cplx(c, v.zEcc);
            c.RFsizemap = max(c.slope.* abs(c.v.ECC) + c.intercept, c.min);

            if ~isfield(c, 'cropPix')
                c.cropPix  = c.v.ANG;
                c.cropPix(c.v.ECC>max([max(v.visfieldHeight) max(v.eccList)]))=NaN;
                if abs(min(c.cortexLength))<abs(max(c.cortexLength))
                    c.cropPix(c.X<0) = NaN;
                elseif abs(min(c.cortexLength))>abs(max(c.cortexLength))
                    c.cropPix(c.X>0) = NaN;
                end
            end
        end

        function v = define_visualmap(v)
            % defines the height, width, resoltion etc of the visual field that is being simulated
            %
            % commented 12/7/2024 IF
            % commented 6/30/2025 ES

            % set defaults
            if ~isfield(v, 'visfieldHeight'); v.visfieldHeight = [-30 30]; end
            if ~isfield(v, 'visfieldWidth'); v.visfieldWidth = [-30 30]; end
            if ~isfield(v,'pixperdeg');     v.pixperdeg = 7;       end
            if ~isfield(v, 'drawthr');     v.drawthr = 1;    end % assumes patients draw percepts when brightness>1

            % generate 1D array 
            v.x = linspace(v.visfieldWidth(1),v.visfieldWidth(2), (v.visfieldWidth(2)-v.visfieldWidth(1)).*v.pixperdeg);
            v.y = linspace(v.visfieldHeight(1),v.visfieldHeight(2), (v.visfieldHeight(2)-v.visfieldHeight(1)).*v.pixperdeg);

            % define 2D grid
            [v.X,v.Y] = meshgrid(v.x, v.y);

            % Make the grid in retinal coordinates
            if ~isfield(v, 'angList');   v.angList = -90:45:90;    end
            if ~isfield(v, 'eccList');  v.eccList = [1 2 3 5 8 13 21 34];    end
            v.gridColor = [1 1 0];
            v.n = 201;
        end

        function [c, v] = define_electrodes(c, v)
            %  Takes in the position of the electrodes (can take multiple) in visual
            %  co-ordinates and pops them onto the cortical surface
            %  Note that there is a 'sister function' c2v_define_electrodes
            %  that takes electrodes on the cortical surface and pops them
            %  onto the visual space
            %
            %  commented 12/7/2024 IF
            %  edited 6/26/2025 ES
            %  commented 6/26/2025 ES

            idx = 1:length(v.e);
            if ~isfield(c, 'e') || ~isfield(c.e, 'radius')
                for ii=1:length(idx);        c.e(idx(ii)).radius = 500/1000;       end
            end
            if ~isfield(c.e, 'shape')
                for ii=1:length(idx);     c.e(idx(ii)).shape = 'round';    end
            end
            if ~isfield(v.e, 'ang')  % if putting in x, y co-ordinates rather than ang and ecc which is the default
                for ii = 1:length(idx)
                    [a, e]= cart2pol(v.e(idx(ii)).x, v.e(idx(ii)).y);
                    v.e(idx(ii)).ang = a*180/pi;
                    v.e(idx(ii)).ecc = e;
                end
            end

            % convert angle and eccentricity to x and y coordinates
            [x_all, y_all] = pol2cart([v.e.ang]*pi/180, [v.e.ecc]);
            for ii = 1:length(idx)
                v.e(ii).x = x_all(ii); v.e(ii).y = y_all(ii);
            end

           % compute area of electrodes
           % convert electrode visual field coordinates to cortical coordinates
           areas = pi * [c.e.radius].^2;
           vx = [v.e.x]; vy = [v.e.y];
           [cx, cy] = p2p_c.v2c_real(c, vx, vy);
           for ii = 1:length(idx)
               c.e(ii).area = areas(ii);
               c.e(ii).x = cx(ii);
               c.e(ii).y = cy(ii);
           end
        end

        function v = c2v_define_electrodes(c,v)
            % If you've defined electrodes (can take multiple) on cortex, this projects them
            % into visual space. Note that there is a 'sister function' define_electrodes
            %  that takes electrodes in visual space and pops them
            %  onto the cortical surface
            %
            %  commented 12/7/2024 IF
            %  edited 6/26/2025 ES
            %  commented 6/26/2025 ES

            %
            idx = 1:length(c.e);

            % map cortical to visual 
            [vx, vy] = p2p_c.c2v_real(c,[c.e.x],[c.e.y]);

            % convert to polar coordinates
            [angs, eccs] = cart2pol(vx, vy);
            angs = angs * 180 / pi;

            % assign back to v.e
            for ii = 1:length(idx)   
                v.e(ii).x = vx(ii);
                v.e(ii).y = vy(ii);
                v.e(ii).ang = angs(ii);
                v.e(ii).ecc = eccs(ii);
            end
        end

        function [v, c] = generate_corticalelectricalresponse(c, v)
            % generates a percept for each electrode by taking the sum of receptive 
            % fields activated by an electrode
            % (weighed by the current at that point on the cortical surface).
            % The percept is normalized so the max is 1
            %
            % commented 12/7/24 IF
            % edited 6/27/25 ES

            % set default model
            if ~isfield(c, 'rfmodel')
                c.rftype = 'ringach';
            end
              
            % generates the sum of weighted receptive fields activated by an electrode
            % normalized so the max is 1
            idx = 1:length(c.e);
            for ii = 1:length(idx) % for each electrode
               
                % get ef and position of current electrode
                ef = double(c.e(ii).ef);
                ex = c.e(ii).x;
                ey = c.e(ii).y;
               
                if min(ex)-1<min(c.X(:)) || max(ex)+1>max(c.X(:)) ... 
                        || min(ey)-1<min(c.Y(:))  ||  max(ey)+1>max(c.Y(:))
                    errordlg('Electrode is either outside or too close to the edge of the cortical sheet');
                else
                    % initialize rf map 
                    rfmap = zeros([size(v.X), 2]);       
                    ctNum = 500000;

                   % identify valid cortical pixels for this electrode     
                    ef_vals = ef(:);
                    crop_vals = c.cropPix(:);
                    valid_mask = ~isnan(crop_vals) & abs(ef_vals) > c.efthr * 255;
                    valid_pix = find(valid_mask);
                    ct = length(valid_pix);

                    if ct <c.e(ii).radius*10
                        if ct ==0;    disp('WARNING! No pixels passed ef threshold.');
                        else;    disp('WARNING! Very few pixels passed ef threshold.');
                        end
                        disp(' try checking the following:');
                        disp('lowering c.efthr or increase stimulation intensity');
                        disp('checking location of electrodes relative visual map');
                        disp('check the sampling resolution of cortex is not too low');
                    end

                    RF_all = p2p_c.generate_corticalcell(ef_vals(valid_pix), valid_pix, c, v);  % [y, x, eye, polarity, nPix]

                    for p = 1:length(valid_pix) % for valid each cortical location:
                       
                        % progress logging
                        %if pixNum>ctNum*2 && (pixNum/ctNum == round(pixNum/ctNum))
                          %  per =round(100*pixNum/length(c.X(:)));
                          %  disp(['generating cortical electrical response ', num2str(per), '% complete']);
                        %end

                        RF = squeeze(RF_all(:, :, :, :, p));  % [y, x, eye, polarity]                

                        % sanity check 
                        if any(isnan(RF(:)))
                            disp('wtf')
                        end
                         
                        % accumulate ON and OFF components
                        rfmap(:, :, 1)  =   rfmap(:, :, 1) + RF(:, :, 1);
                        rfmap(:, :, 2)  =   rfmap(:, :, 2) + RF(:, :, 2);  
                    end

                    % normalize
                     v.e(idx(ii)).rfmap = rfmap./max(abs(rfmap(:)));
                end                        
            end
        end


        % function [v, c] = generate_corticalvisualresponse_spatial(c, v)
        %     if ~isfield(c, 'rfmodel')
        %         c.rftype = 'ringach';
        %     end
        % 
        %     c.target.img = NaN(size(c.X(:)));
        %     % convolve the of  receptive fields with the visual stimulus
        %     for pixNum = 1:length(c.X(:)) % for each cortical location
        %         if (pixNum/100000 ==round(pixNum/100000))
        %             per = round(100*pixNum/length(c.X(:)));
        %             disp(['generating cortical visual response ', num2str(per), '% complete']);
        %         end
        %         if ~isnan(c.cropPix(pixNum))
        %             RF = p2p_c.generate_corticalcell(1, pixNum, c, v);
        %             if ~isempty(find(~isnan(RF(:)), 1))
        %                 on =  squeeze(RF(:, :, 1, 1)).*v.target.img; % on response to target
        %                 off =  squeeze(RF(:, :, 1, 2)).*v.target.img; % suppressive response to target
        %                 c.target.img(pixNum) =sum(on(:))-sum(off(:));
        %             else
        %                 c.target.img(pixNum) = 0;
        %             end
        %         end
        %     end
        %     c.target.img = reshape(c.target.img, size(c.X)); % the pattern of activity in space
        % end

        % function [c, v] = generate_visualresponse_temporal(tp, c, v)
        %     % leaky integrator of stimulus attracted to a baseline which shifts over
        %     % time.
        % 
        %     t = 0:tp.dt:(v.target.dur-tp.dt);
        %     R(1) = 0;
        %     v.target.tc = zeros(size(t));
        %     v.target.tc(t>=v.target.on & t<= v.target.off) = v.target.amp;
        %     delayFrames = round(tp.vis.delay./tp.dt);
        %     c.target.response = zeros(size(t));  b = zeros(size(t));
        %     dy2 = 0; dy = 0; db = 0;
        %     for i=(delayFrames+1):(length(t)-1)
        %         if v.target.tc(i-delayFrames) > 0  % stimulus is on - adapt
        %             db = (v.target.tc(i-delayFrames)-b(i))/tp.vis.tau_b_on;
        %         else  % stimulus is off - recover
        %             db = (v.target.tc(i-delayFrames)-b(i))/tp.vis.tau_b_off;
        %         end
        %         dy = (v.target.tc(i-delayFrames) - c.target.response(i)-b(i))/tp.vis.tau_1;
        %         b(i+1) = b(i)+db*tp.dt;
        %         c.target.response(i+1) = c.target.response(i)+dy*tp.dt;
        %     end
        %     v.target.t = t;
        % end

        % function G = Gauss_2D(v,x0,y0,theta,sigma_x,sigma_y)
        %     Generates oriented 2D Gaussian on meshgrid v.X,v.Y
        %     aa = cos(theta)^2/(2*sigma_x^2) + sin(theta)^2/(2*sigma_y^2);
        %     bb = -sin(2*theta)/(4*sigma_x^2) + sin(2*theta)/(4*sigma_y^2);
        %     cc = sin(theta)^2/(2*sigma_x^2) + cos(theta)^2/(2*sigma_y^2);
        %     G = exp( - (aa*(v.X-x0).^2 + 2*bb*(v.X-x0).*(v.Y-y0) + cc*(v.Y-y0).^2));
        % end

        function RF = generate_corticalcell(ef, pix_idx, c, v)
            try
                % Get pixel-wise parameters
                x0 = c.v.X(pix_idx);      % [nPix x 1]
                y0 = c.v.Y(pix_idx);      % [nPix x 1]
                od = c.ODmap(pix_idx);    % [nPix x 1]
                theta = pi - c.ORmap(pix_idx);  % [nPix x 1]
                sigma_x = c.RFsizemap(pix_idx) * c.ar;
                sigma_y = c.RFsizemap(pix_idx);
        
                X = v.X; Y = v.Y;
                [nY, nX] = size(X);
                nPix = numel(pix_idx);
        
                RF = zeros(nY, nX, 2, 2, nPix);  % Preallocate for Ringach
                model = c.rfmodel;
        
                % Repeat coordinate grid
                Xr = reshape(X, 1, 1, []); Yr = reshape(Y, 1, 1, []);
        
                for p = 1:nPix
                    dx = X - x0(p); dy = Y - y0(p);
        
                    aa = cos(theta(p))^2 / (2*sigma_x(p)^2) + sin(theta(p))^2 / (2*sigma_y(p)^2);
                    bb = -sin(2*theta(p)) / (4*sigma_x(p)^2) + sin(2*theta(p)) / (4*sigma_y(p)^2);
                    cc = sin(theta(p))^2 / (2*sigma_x(p)^2) + cos(theta(p))^2 / (2*sigma_y(p)^2);
        
                    switch model
                        case 'scoreboard'
                            G = ef(p) * exp(-((dx.^2)/0.0001 + (dy.^2)/0.00001));
                            RF(:, :, 1, 1, p) = G;
                            RF(:, :, 2, 1, p) = G;
        
                        case 'smirnakis'
                            G = ef(p) * exp(-(aa*dx.^2 + 2*bb*dx.*dy + cc*dy.^2));
                            G = G ./ sum(G(:));
                            RF(:, :, 1, 1, p) = od(p) * G;
                            RF(:, :, 2, 1, p) = (1 - od(p)) * G;
        
                        case 'ringach'
                            tmp = exp(-(aa*dx.^2 + 2*bb*dx.*dy + cc*dy.^2));
                            A = sqrt(sum(tmp(:) > 0.2) / v.pixperdeg^2);
                            d = c.DISTmap(pix_idx(p)) * A;
        
                            wplus = c.ONOFFmap(pix_idx(p));
                            wminus = 1 - wplus;
        
                            x_on = x0(p) - (d / 2) * cos(theta(p));
                            y_on = y0(p) + (d / 2) * sin(theta(p));
                            x_off = x0(p) + (d / 2) * cos(theta(p));
                            y_off = y0(p) - (d / 2) * sin(theta(p));
        
                            hplus_on = exp(-(aa*(X - x_on).^2 + 2*bb*(X - x_on).*(Y - y_on) + cc*(Y - y_on).^2));
                            hplus_off = exp(-(aa*(X - x_off).^2 + 2*bb*(X - x_off).*(Y - y_off) + cc*(Y - y_off).^2));
        
                            hplus_on = (hplus_on ./ max(abs(hplus_on(:)))) * wplus;
                            hplus_off = (hplus_off ./ abs(max(hplus_off(:)))) * wminus * c.onoff_ratio;
        
                            hminus_off = -0.4 * hplus_off;
                            hminus_on = -0.4 * hplus_on;
        
                            excitatory = hplus_on - c.onoff_ratio * hplus_off;
                            inhibitory = hminus_on - c.onoff_ratio * hminus_off;
        
                            RF(:, :, 1, 1, p) = od(p) * ef(p) * excitatory;
                            RF(:, :, 2, 1, p) = (1 - od(p)) * ef(p) * excitatory;
                            RF(:, :, 1, 2, p) = od(p) * ef(p) * inhibitory;
                            RF(:, :, 2, 2, p) = (1 - od(p)) * ef(p) * inhibitory;
                    end
                end
        
            catch ME
                warning('Falling back to per-pixel loop due to memory constraints or error: %s', ME.message);
                RF = zeros(size(v.X,1), size(v.X,2), 2, 2, numel(pix_idx));
                for k = 1:numel(pix_idx)
                    RF(:, :, :, :, k) = p2p_c.generate_corticalcell(ef(k), pix_idx(k), c, v);
                end
            end
        end

        function [trl_array,v] = generate_phosphene(v, tp, trl)
            % finds the phosphene trl.max_phosphene corresponding to the brightest moment in
            % time and estimates trl.ellipse which represents the patient
            % drawings 
            %
            % arguments : 
            % v - struct of visual space information
            % tp - temporal parameters
            % trl - single trial parameters struct
            %
            % return values: 
            % trl_array - a struct array of length # electrodes that now
            % contains phosphene info corresponding to the brightest moment
            % and estimates trl.ellipse which represents the patient
            % drawings for each electrode
            % v - unchanged
            %
            % commented 12/7/24 IF
            % edited 7/24/24 ES


            % Ensure all fields that will be added exist in trl
            required_fields = {'ptid', 'spikestrength', 'max_phosphene', ...
                   'sim_area', 'ellipse', 'sim_brightness'};
     
          
            for f = 1:length(required_fields)
                field = required_fields{f};
                if ~isfield(trl, field)
                     trl.(field) = [];  % initialize with empty
                end
            end

            
            n_electrodes = length(v.e); % number of electrodes passed
            trl_array = repmat(trl, 1, n_electrodes); % preallocate output array 

            for e = 1:n_electrodes

                curr_trl = trl;
                curr_trl.e = e;
                
                if isnan(trl.freq)
                    curr_trl.max_phosphene = v.e(curr_trl.e).rfmap; 
                    curr_trl.spikestrength = 1; % the scaling due to current integration
                else
                    curr_trl = p2p_c.convolve_model(tp, curr_trl);
                    curr_trl.max_phosphene = v.e(curr_trl.e).rfmap.*max(curr_trl.spikestrength);
                end

                % calculate the size of the image
                curr_trl.sim_area = (1/v.pixperdeg.^2) * sum(curr_trl.max_phosphene(:) > v.drawthr)/2; % calculated area of phosphene based on mean of left and right eyes
              
                if ~isempty(curr_trl.max_phosphene)
                    for i=1:2 % left and right eye
                        p = p2p_c.fit_ellipse_to_phosphene(curr_trl.max_phosphene(:,:,i)>v.drawthr,v);
                        curr_trl.ellipse(i).x = p.x0;
                        curr_trl.ellipse(i).y = p.y0;
                        curr_trl.ellipse(i).sigma_x = p.sigma_x;
                        curr_trl.ellipse(i).sigma_y = p.sigma_y;
                        curr_trl.ellipse(i).theta = p.theta;
                    end
                    % what rule to use to translate phosphene image to brightness?
                    beta = 6; % soft-max rule across pixels for both eyes
                    curr_trl.sim_brightness = ((1/v.pixperdeg.^2) * sum(curr_trl.max_phosphene(:).^beta)^(1/beta));  % IF CHECK
                else
                    curr_trl.sim_brightness = [];
                end
                trl_array(e) = curr_trl;
            end
        end

        function [tp, nc] = define_temporalparameters(varargin)
            % defines the parameters of the temporal model, defaults are
            % based on Fine and Boynton, 2024
            % https://www.nature.com/articles/s41598-024-65337-1#citeas
            %
            % commented 12/7/2024 IF

            if nargin<1
                tp = [];
            else
                tp = varargin{1};
            end

            % parameters controlling the first stage of rapid integration
            % of current by cell, describes the effect of changing pulse
            % width
            if ~isfield(tp, 'dt');       tp.dt = .001 * 10^-3; end % Represents time sampling in ms, should be no larger than 1/10 of tau1
            if ~isfield(tp, 'tau1');   tp.tau1 = 0.0003; end    % First stage temporal integrator fixed based on Nowak and Bullier, 1998, 10.1007/s002210050304
            if ~isfield(tp, 'refrac');   tp.refrac = 100; end  % Extent of attenuation for spiking refractory period
            if ~isfield(tp, 'delta');   tp.delta = 0.001; end  % Decay parameter for refractory period

            % parameters controlling slower second stages of integration,
            % describes the effect of changing the frequency, this model
            % isn't great
            if ~isfield(tp, 'tSamp');   tp.tSamp = 1000; end % Subsampling to speed things up
            if ~isfield(tp, 'tau2');   tp.tau2 =  0.025; end  % Slower second stage of integration
            if ~isfield(tp, 'ncascades');  tp.ncascades = 3;   end % number of cascades in the slow filter of the temporal convolution
            if ~isfield(tp, 'gammaflag');   tp.gammaflag = 1;   end            %  include second stage gamma


            % leak out of charge accumulation
            %   tp.flag_cl=0; % 1 if you want to charge to leak back out of the system
            %   tp.tau2_cl = tp.tau2_ca * 8; % used for the conv model, fe uses p.tau2_ca

            % nonlinearity response parameters
            if ~isfield(tp, 'spikemodel')
                tp.spikemodel = 'convolve';
            end

            if strcmp(tp.spikemodel, 'convolve')
                if ~isfield(tp, 'model');   tp.model = 'compression';   end
                if ~isfield(tp, 'sc_in');   tp.sc_in = 0.5663;   end
            elseif strcmp(tp.spikemodel, 'integratefire')
                if ~isfield(tp, 'model');   tp.model = 'compression';   end
                if ~isfield(tp, 'sc_in');   tp.sc_in = .5663;  end % 10000/3;   end
            else
                error('model variant not defined, should be "compression" or "linear"');
            end

            if strcmp(tp.model, 'compression')
                if ~isfield(tp, 'power'); tp.power =  15.5901; end % chosen cos max brightness rating
                if ~isfield(tp, 'sc_out'); tp.sc_out = 10; end % fit using Winawer brightness data
            elseif strcmp(tp.model, 'sigmoid')
                disp('using sigmoid semisaturation constant')
                tp.asymptote = 2000;
                tp.e50 = 500; % electrical semisaturation constant
            elseif strcmp(tp.model, 'normcdf')
                disp('using normcdf semisaturation constant')
                tp.asymptote = 1500;
                tp.mean = 750;
                tp.sigma = 175;
            elseif strcmp(tp.model, 'weibull')
                disp('using weibull semisaturation constant')
                tp.asymptote = 1000;
                tp.thresh = 600;
                tp.beta = 3.5;
            end
        end

        %% psychophysics
        % function v = generate_visualtarget(v)
        %     [x, y] = pol2cart(v.target.ang*pi/180, v.target.ecc);
        %     v.target.img = sqrt(((v.X-x).^2)+((v.Y-y).^2))<v.target.rad;
        % end

        % function [err, thresh] = loopall_find_threshold(tp,T)
        %     %
        %     % Runs the 'conv' model to get thresholds based on trials in the table 'T'.
        %     % returns the SSE and thresholds.  Table must contain fields holding the
        %     % following parameters for each trial:
        %     % pw      pulse width (sec)
        %     % dur     trial duration (sec)
        %     % freq    pulse frequency (Hz)
        %     % amp     amplitude at detection threshold
        %     if ~isfield(tp, 'nReps')
        %         tp.nReps = 12;
        %     end
        %     thresh = NaN(1,size(T,1));
        %     for i=1:size(T,1)
        %         % define trial parameters based on values in the table
        %         clear trl;  trl.pw = T.pw(i);   trl.amp = 1;    trl.dur = T.dur(i);     trl.freq = T.freq(i);   trl.simdur = 3; %sec
        %         trl = p2p_c.define_trial(tp,trl);
        % 
        %         tp.sc_in_ind = 1;     % separate 'scFac' for each experiment
        %         if isfield(tp,'experimentList')  % set tp_thresh accordingly for this trial
        %             experimentNum = find(strcmp(T.experiment{i},tp.experimentList));
        %             if ~isempty(experimentNum)  % set thresh_resp for this experiment.
        %                 tp.sc_in_ind = experimentNum;
        %             end
        %         end
        %         thresh(i)= p2p_c.find_threshold(trl,tp);
        %     end
        % 
        %     err = nansum((thresh-T.amp').^2);
        %     %    disp(sprintf('tau1 = %g, tau2 = %g, power = %5.2f err= %5.4f  sc_in= %5.4f',  tp.tau1,tp.tau2,tp.power,err, tp.sc_in));
        %     disp(fprintf('mean = %g, sigma = %g,  err= %5.4f\n',  tp.mean,tp.sigma,err));
        %     if isfield(tp,'experimentList')
        %         for i = 1:length(tp.experimentList)
        %             disp(fprintf('%10s: %g\n',tp.experimentList{i},tp.sc_in(i)));
        %         end
        %     end
        % end
        % 
        % function [err, thresh] = loop_find_threshold(tp,T, varargin)
        % 
        %     %
        %     % Runs the 'conv' model to get thresholds based on trials in the table 'T'.
        %     % returns the SSE and thresholds.  Table must contain fields holding the
        %     % following parameters for each trial:
        %     % pw      pulse width (sec)
        %     % dur     trial duration (sec)
        %     % freq    pulse frequency (Hz)
        %     % amp     amplitude at detection threshold
        %     if ~isfield(tp, 'nReps')
        %         tp.nReps = 12;
        %     end
        %     thresh = NaN(size(T,1), 1);
        %     for i=1:size(T,1)
        %         % define trial parameters based on values in the table
        %         clear trl;  trl.pw = T.pw(i);   trl.amp = 1;    trl.dur = T.dur(i);     trl.freq = T.freq(i);   trl.simdur = 3; %sec
        %         trl = p2p_c.define_trial(tp,trl);
        % 
        %         if nargin==2
        %             thresh(i)= p2p_c.find_threshold(trl,tp);
        %         elseif nargin ==3
        %             thresh(i)= p2p_c.find_threshold(trl,tp, varargin{1});
        %         end
        %     end
        % 
        %     if strcmp('amp',T.Properties.VariableNames)
        %         err = nansum((thresh-T.amp').^2);
        %         disp(['err = ',  num2str(round(err, 3))]);
        %     else
        %         err = NaN;
        %     end
        % 
        %     if isfield(tp,'experimentList')
        %         for i = 1:length(tp.experimentList)
        %             disp(sprintf('%10s: %g',tp.experimentList{i},tp.sc_in(i)));
        %         end
        %     end
        % end

        % function err = fit_brightness(tp, T, varargin)
        %     if nargin==3
        %         nc = varargin{1};
        %         [loop_trl] = p2p_c.loop_model(tp,T, nc);
        %     else
        %         [loop_trl] = p2p_c.loop_model(tp,T);
        %     end
        % 
        % 
        %     y_est = [loop_trl.max_temporal_response]; y = [T.brightness];
        %     y_est = reshape(y_est, length(y), 1);
        %     y = reshape(y, length(y), 1);
        %     ind = ~isnan(y_est) & ~isnan(y);
        %     err = sum((y(ind)- y_est(ind)).^2);
        %     disp(sprintf('tau1 =%5.4f, tau2 =%5.4f, power =%5.4f, sc_in =%5.4f, sc_out =%5.4f, sse = %5.4f',  ...
        %         tp.tau1, tp.tau2, tp.power,tp.sc_in, tp.sc_out, err));
        % end

        % function [loop_trl] = loop_model(tp,T, varargin)
        %     %
        %     % Runs the 'conv' model to get thresholds based on trials in the table 'T'.
        %     % returns the SSE and thresholds.  Table must contain fields holding the
        %     % following parameters for each trial:
        %     % pw      pulse width (sec)
        %     % dur     trial duration (sec)
        %     % freq    pulse frequency (Hz)
        %     % amp     amplitude at detection threshold
        %     if nargin==3
        %         nc = varargin{1};
        %     end
        %     for i=1:size(T,1)
        %         % define trial parameters based on values in the table
        %         clear trl;  trl.pw = T.pw(i);   trl.amp = T.amp(i);    trl.dur = T.dur(i);     trl.freq = T.freq(i);   trl.simdur = 3; %sec
        %         trl = p2p_c.define_trial(tp,trl);
        % 
        %         % define impulse response
        %         if isfield(tp,'tSamp')
        %             if tp.tSamp~=1% down-sample the time-vectors
        %                 t = trl.t(1:tp.tSamp:end);
        %             end
        %         else
        %             t = trl.t;
        %         end
        %         dt = t(2)-t(1);
        %         h = p2p_c.gamma(tp.ncascades,tp.tau2,t);            % Generate the n-cascade impulse response
        %         tid = find(cumsum(h)*dt>.999,1,'first'); % Shorten the filter if needed to speed up the code.
        %         if ~isempty(tid)
        %             h = h(1:tid);
        %         else
        %             sprintf('Warning: gamma hdr might not have a long enough time vector');
        %         end
        %         trl.imp_resp = h;  % close enough to use h
        %         if strcmp(tp.spikemodel, 'integratefire')
        %             if exist('nc')
        %                 loop_trl(i) = p2p_c.spike_model(tp, trl, nc);
        %             else
        %                 loop_trl(i) = p2p_c.spike_model(tp, trl);
        %             end
        %         else
        %             loop_trl(i) = p2p_c.spike_model(tp, trl);
        %         end
        %     end
        % end

        % function [trl, nc,v] = integratefire_model(tp, trl, nc)
        % 
        %     % Implements the 'Adaptive Exponential Integrate-and-fire' model (AdEx)
        %     % which seems to be the gold standard for a flexible model for neural
        %     % spiking.  Lots of references point to pretty much the same model,
        %     % including the one Ione found:
        %     %
        %     % https://neuronaldynamics.epfl.ch/online/Ch6.S1.html
        %     %
        %     % I've also added in the 'absoulte refractory period' (parameter tARP).
        %     %
        %     % This simulates the response to a step function with amplitude stimAmp and
        %     % duration stimDur.
        %     %
        %     % By messing with parameters, you can apparently get the model to do
        %     % neuronally thinks like tonic, bursting, delays, transience, and
        %     % adaptation.  Here's a video of a humorless guy describing the behaviors
        %     % (but not which parameters do what).  Something about 'phase-plane'
        %     % analysis can help.
        %     %
        %     % https://www.youtube.com/watch?v=u8qwboY-zA4
        %     %
        %     % https://www.youtube.com/watch?v=GbXAf0iMvzs
        % 
        %     % stimulation environment parameters
        % 
        % 
        % 
        %     nt = length(trl.t);
        %     dt = tp.dt;
        %     Iinput = trl.pt; % *1000000;  % picoamps to microamps
        % 
        %     % zero out vector to store membrane potential (v) and adapation variable (w)
        %     for cc = 1:length(nc)
        %         clear names
        %         % Hack:  turns out accessing values in a structure is about 10x slower than
        %         % accessing a 'loose' variable.  Who knew?  So this pulls all fields out
        %         % of a structure for use in the loop.
        % 
        %         names = fieldnames(nc(cc)); % take everything out of the structure to speed things up
        %         for i=1:length(names)
        %             eval(sprintf('%s = nc(%d).%s;', names{i},cc,names{i}))
        %         end
        %         nARP = round(tARP /dt); % number of time points where there is absolute refactory period
        % 
        %         d = 0;
        %         v = zeros(nt,1); v(1) = El;
        %         w = zeros(nt,1);
        %         for timei = 1:(nt-1)  % simulation in loop over time
        %             v(timei+1) = v(timei) + dt*(-gl*(v(timei)-El) + ...
        %                 gl*delT*exp((v(timei)-vt)/delT) + ...
        %                 Iinput(timei)*1000 - w(timei) )*R;
        % 
        %             % v(timei+1) = v(timei) + dt*(-gl*(v(timei)-El) +Iinput(timei))*R;
        %             if d>0
        %                 v(timei+1) = vreset;
        %             end
        %             % note: updat e w by pre-updated v
        %             w(timei+1) = w(timei) + dt*(a*(v(timei)-El)-w(timei))/tauW;
        % 
        %             % check for a spike
        %             if v(timei+1) >= vpeak
        %                 v(timei) = vpeak;
        %                 v(timei+1) = vreset;
        %                 w(timei+1)  = w(timei+1)  + b;
        %                 d = nARP;
        %             end
        %             d = max(d-1,0);
        %         end
        %         indexFire(:, cc) = v == nc(cc).vpeak; % a matrix showing when each cell fired
        %     end
        %     if length(nc)>1
        %         trl.ptid = find(sum(indexFire,2));
        %         trl.spikestrength =  sum(indexFire(trl.ptid,:),2)./length(nc);
        %     else
        %         trl.ptid = find(indexFire);
        %         trl.spikestrength = indexFire(trl.ptid);
        %     end
        % end

        % function trl = membrane_depolarization(tp, trl, nc)
        %     % models neural suppression as a result of electrical stimulation
        % 
        %     if ~isfield(nc, 'tauMDP'); nc.tauMDP = 0.25; end % time course of suppression
        %     if ~isfield(nc, 'strengthMDP'); nc.strengthMDP = 1; end % strength of suppression
        % 
        %     Iext = nc.Vext.*tp.dt; % input is scaled by Ifac and resistance
        % 
        %     trl.MDP = zeros(size(nc.Vext)); % set initial membrane depolarization to 0
        %     Iext(Iext<0) = 0; % rectification, here we assume suppression based on cathodic component, but likely anodic?
        % 
        %     for c = 1:length(trl.pt) - 1
        %         % calculate membrane depolarization over time using a leaky
        %         % integrator
        %         trl.MDP(:, c+1) =  (trl.MDP(:, c) + Iext(:, c)) - ((tp.dt/nc.tauMDP).*trl.MDP(:, c)); % assume suppression driven by charge
        %     end
        %     trl.MDP = trl.MDP.*nc.strengthMDP;
        % end

        % function threshold = find_threshold(trl, tp, varargin)
        %     % Find amplitudes at threshold with the leaky integrate and
        %     % fire model
        %     % takes in trial, tp, and optional fitParams
        %     % finds and returns the trl.amp for which the max output of the
        %     % model for that trial, trial.resp, is equal to fitParams.thr
        % 
        %     if nargin ==2
        %         nCells = 1;
        %         nc = p2p_c.define_integratefirecells(nCells);
        %     elseif nargin == 3 && ~isempty(varargin{1})
        %         if isstruct(varargin{1})
        %             nc = varargin{1};
        %         elseif strcmp(tp.spikemodel, 'integratefire')
        %             nCells = varargin{1};
        %             nc = p2p_c.define_integratefirecells(nCells);
        %         end
        %     end
        %     if nargin == 4
        %         verbose = varargin{2};
        %     else
        %         verbose = 0;
        %     end
        %     if ~isfield(tp, 'nReps')
        %         tp.nReps = 12;
        %     end
        % 
        %     % first find the lowest 'hi' response
        % 
        %     resp = 0;
        %     if verbose;  disp('finding upper amplitude for search'); end
        %     trl.amp = .5;
        % 
        %     while resp<tp.thresh_resp
        %         trl.amp = trl.amp*2;
        %         trl = p2p_c.define_trial(tp,trl);
        %         if strcmp(tp.spikemodel, 'integratefire')
        %             trl = p2p_c.spike_model(tp, trl, nc);
        %         elseif strcmp(tp.spikemodel, 'convolve')
        %             trl = p2p_c.spike_model(tp, trl);
        %         else
        %             error('spiking model not recognized');
        %         end
        %         resp = max(trl.resp);
        %         disp(['amp = ', num2str(trl.amp), ' resp = ', num2str(resp), ' thresh = ', num2str(tp.thresh_resp)] )
        %     end
        % 
        %     % 3/7/25 gmb changed lo from 0 to trl.amp/2 since we know thresh is between amp/2 and amp
        %     hi = trl.amp;  
        %     if trl.amp > 1
        %         lo = trl.amp/2;
        %     else
        %         lo = 0;
        %     end
        %     % then do the binary search
        %     for i = 1:tp.nReps
        %         trl.amp  = (hi+lo)/2;
        %         trl = p2p_c.define_trial(tp,trl);
        %         if strcmp(tp.spikemodel, 'integratefire')
        %             trl = p2p_c.spike_model(tp, trl, nc);
        %         elseif strcmp(tp.spikemodel, 'convolve')
        %             trl = p2p_c.spike_model(tp, trl);
        %         else
        %             error('spiking model not recognized');
        %         end
        %         resp = max(trl.resp);
        %         if max(resp(:)) > tp.thresh_resp
        %             hi = trl.amp;
        %         else
        %             lo = trl.amp;
        %         end
        %         disp(['amp = ', num2str(trl.amp), ' resp = ', num2str(resp), ' thresh = ', num2str(tp.thresh_resp)] )
        %     end
        %     threshold = trl.amp;
        % end


        % function nc = define_integratefirecells(varargin)
        % 
        %     if nargin == 1 & isstruct(varargin{1})
        %         nc = varargin{1};
        %         nCells = length(nc);
        %         for cc = 1:nCells
        %             if ~isfield(nc(cc), 'gl') || isempty(nc(cc).gl);      nc(cc).gl =  10; end % leak conductance (10)
        %             if ~isfield(nc(cc), 'El') || isempty(nc(cc).El);      nc(cc).El = -58; end % leak voltage (-58)
        %             if ~isfield(nc(cc), 'vt') || isempty(nc(cc).vt);      nc(cc).vt  = -50; end % resting membrane potential (-50)
        %             if ~isfield(nc(cc), 'delT') || isempty(nc(cc).delT);  nc(cc).delT =   2; end % spike Width Factor (2)
        %             if ~isfield(nc(cc), 'a') || isempty(nc(cc).a);        nc(cc).a  =   2; end % resonator/integrator variable (2)
        %             if ~isfield(nc(cc), 'tauW') || isempty(nc(cc).tauW);  nc(cc).tauW = 0.15; end % adaptation decay time (0.12)
        %             if ~isfield(nc(cc), 'b') || isempty(nc(cc).b);        nc(cc).b = 100; end % adaptation jump (update for w) (100)
        %             if ~isfield(nc(cc), 'vreset') || isempty(nc(cc).vreset); nc(cc).vreset = -46; end % voltage reset (-46)
        %             if ~isfield(nc(cc), 'vpeak') || isempty(nc(cc).vpeak);    nc(cc).vpeak = 0; end % spike threshold
        %             if ~isfield(nc(cc), 'tARP') || isempty(nc(cc).tARP);      nc(cc).tARP = 0 ; end % 5*10^-3;  % absolute refractory period
        %         end
        %     else
        %         if nargin==0
        %             nCells = 1;
        %         elseif nargin == 1 && isnumeric(varargin{1})
        %             nCells = varargin{1}; % using the number of cells and filling in values.
        %         end
        % 
        %         for cc = 1:nCells  % populate any missing values
        %             if nCells == 1
        %                 nc(1).R = 5;
        %             end
        %             nc(cc).R = 4.5+rand(1);
        %             nc(cc).gl =  10;  % leak conductance (10)
        %             nc(cc).El = -58;  % leak voltage (-58)
        %             nc(cc).vt  = -50;  % resting membrane potential (-50)
        %             nc(cc).delT = 2;  % spike Width Factor (2)
        %             nc(cc).a  = 2;  % resonator/integrator variable (2)
        %             nc(cc).tauW = 0.15;  % adaptation decay time (0.12)
        %             nc(cc).b = 100;  % adaptation jump (update for w) (100)
        %             nc(cc).vreset = -46;  % voltage reset (-46)
        %             nc(cc).vpeak = 0;  % spike threshold
        %             nc(cc).tARP = 0 ; % 5*10^-3;  % absolute refractory period
        % 
        %         end
        %     end
        % end


        % function [trl, nc] = spike_model(tp, trl, varargin)
        %     % There are two variations for predicting neural / perceptual responses to
        %     % electrical stimulation
        %     % convolve and integrate fire
        %     % takes in tp, trl and variable argument nc if using an
        %     % integrate and fire model
        % 
        %     if ~isfield(tp, 'spikemodel')
        %         tp.spikemodel = 'convolve';
        %     end
        %     if strcmp(tp.spikemodel, 'convolve')
        %         trl = p2p_c.convolve_model(tp, trl);
        %     elseif strcmp(tp.spikemodel, 'integratefire')
        %         if nargin>2
        %             nc = p2p_c.define_integratefirecells(varargin{1}); % set up the cells for the integrate fire model
        %         else
        %             nc = p2p_c.define_integratefirecells(1);
        %         end
        %         [trl, nc] = p2p_c.integratefire_model(tp, trl, nc);
        %     end
        %     if tp.gammaflag
        %         trl = p2p_c.slowgamma(tp, trl);
        %     else
        %         trl.resp = zeros(size(trl.pt));
        %         trl.resp(trl.ptid) = trl.spikestrength;
        %         trl.max_temporal_response = max(trl.resp);
        %     end
        % end

        function trl = convolve_model(tp, trl)
            % Implements 'finite_element' using the closed-form solution to
            % the response to a pulse   Can be faster than 'finite_element'.
            % Assumes square pulse trains.
            %
            % tSamp is the temporal sub-sampling factor. Since tau2 is
            % relatively long, we can get away with a coarser temporal
            % sampling for the last convolution stage.  tSamp of 1000 works
            % well. Advise comparing to tSamp = 1 to check for innacuracy.
            % Also advise comparing 'convolve_model' to 'finite_element'
            % model which should be consiered the ground truth.
            %
            % Note: model only returns 'R3', 'spike' and 'resp' as output.
            % R1 and R2 (rectified R1) timecourses are not generated, so
            % the R2 of 'simpleleakyintegrator' has to obtained through the
            % 'finite_element' function.
            %
            % 'tt' is also returned, which is the temporally subsampled 't'
            % vector. Good for plotting 'spike' and 'resp'.

            % written GMB 6/17/2022

            % Assume the pulse train, pt, is a sequence of discrete jumps
            % in current. Find the 'events' where the pulse train, pt,
            % jumps up or down.
            %     Rconvtmp =  zeros(1,tp.ncascades+1); CHECK WITH GEOFF

            % Since spikes are sparse, manually convolve the 'spikes' with
            % the impulse response function at lower temporal
            % resolution

            if isfield(tp,'tSamp') * tp.tSamp~=1 % down-sample the time-vectors
                t = trl.t(1:tp.tSamp:end);
            else
                t = trl.t;
            end
            trl.ptid = find(diff(trl.pt))+1;

            % R will hold the values of R1 at the event times.
            Rtmp = zeros(1,length(trl.ptid));
            wasRising = 0;
            % Loop through the events, calculating R1 at the end of the event
            % and add impulse responses when R1 peaks and is after the refractory period.
            spikeId = logical(size(Rtmp));
            for i=1:(length(trl.ptid)-1)
                %    tNow = trl.t(ptid(i+1));
                delta = trl.t(trl.ptid(i+1))-trl.t(trl.ptid(i));  % time since last 'event'
                % Closed form solution to leaky integrator that predicts
                % R(i+1) from R(i), delta and tau1:
                Rtmp(i+1) = trl.pt(trl.ptid(i))*tp.tau1*(1-exp(-delta/tp.tau1)) + ...
                    Rtmp(i)*exp(-delta/tp.tau1);

                % Add a spike if:
                % (1) R1 is going down since last event
                % (2) R1 was going up before that, and
                % (3) we're past the refractory period since the last spike
                if Rtmp(i+1)<Rtmp(i) && wasRising
                    spikeId(i) = 1; % check spike id identical in both loops IF CHECK
                    wasRising = 0;  % no longer rising
                    %    lastSpikeTime = trl.t(ptid(i+1));
                else
                    wasRising =1;
                end
            end

            R1= Rtmp*1000;
            R1(R1<0)=0;
            if isempty(R1) % no spikes at all
                trl.resp = zeros(1, length(t));
                trl.resp_lin = trl.resp;
                trl.tt = t;  % for plotting
                trl.max_temporal_response = 0;
                trl.spikeWhen= NaN;
                trl.imp_resp = NaN;
                trl.spikestrength = NaN; trl.spikes_norefrac = NaN;
            else
                spikeId (R1<=0) = 0; % sometimes non-spikes identified as spikes
                % pull out only spike events
                trl.spikestrength = R1(spikeId);
                trl.ptid = trl.ptid(spikeId);
            end
        end

        % function trl = slowgamma(tp, trl)
        %     % takes the estimate of spiking (spiking strength for the
        %     % convolve model, spikes across multiple neurons for the
        %     % integrate fire model and passes through a slow filter
        %     % that represents the brain
        % 
        % 
        %     h = p2p_c.gamma(tp.ncascades,tp.tau2,trl.t);            % Generate the n-cascade impulse response
        %     tid = find(cumsum(h)*tp.dt>=.999,1,'first');      % Shorten the filter if needed to speed up the code.
        %     if ~isempty(tid)
        %         h = h(1:tid);
        %         trl.imp_resp = h;  % close enough to use h
        %     else
        %         disp('Warning: gamma hdr might not have a long enough time vector');
        %     end
        %     impFrames = 0:(length(trl.imp_resp)-1);
        %     resp = zeros(1,length(trl.t)+length(trl.imp_resp));        % zero stuff out
        % 
        %     if strcmp(tp.spikemodel, 'convolve')
        %         % add an artificial refractory period, not needed for the
        %         % integrate and fire model
        %         interspike = [1,diff(trl.t(trl.ptid))];
        %         trl.spikes_norefrac = trl.spikestrength;
        %         trl.spikestrength = trl.spikestrength.*(1-exp(-tp.refrac*(interspike+tp.delta)));
        %     end
        % 
        %     for i=1:length(trl.spikestrength)
        %         id = find(trl.t>trl.t(trl.ptid(i)),1,'first');
        %         resp(id+impFrames)  =   ...
        %             resp(id+impFrames) + trl.spikestrength(i)*trl.imp_resp;
        %     end
        % 
        %     trl.resp_lin = resp;
        %     resp = p2p_c.nonlinearity(tp, resp);
        % 
        %     % save the time-course of the response for output
        %     trl.resp = resp(1:length(trl.t));
        %     trl.spikeWhen = trl.ptid;
        %     trl.max_temporal_response = max(trl.resp);
        % end

        %% utilities
        % function out = chronaxie(p,pw)
        %     out = p.amp./(p.tau*(1-exp(-pw/p.tau)));
        % end

        % function y = gamma(n,k,t)
        %     %   y=gamma(n,k,t)
        %     %   returns a gamma function on vector t
        %     %   y=(t/k).^(n-1).*exp(-t/k)/(k*factorial(n-1));
        %     %   which is the result of an n stage leaky integrator.
        % 
        %     %   6/27/95 Written by G.M. Boynton at Stanford University
        %     %   4/19/09 Simplified it for Psychology 448/538 at U.W.
        %     %
        %     y = (t/k).^(n-1).*exp(-t/k)/(k*factorial(n-1));
        %     y(t<0) = 0;
        % end

        % function y = nonlinearity(tp,x)
        %     if ~isfield(tp, 'sc_in');    tp.sc_in = 1;
        %     else
        %         sc_in  =  tp.sc_in;
        %     end
        %     % some of our favorite static nonlinearities:
        %     switch tp.model
        %         case 'sigmoid'
        %             y = sc_in .* x.^tp.power./(x.^tp.power + tp.sigma.^2);
        %         case 'normcdf'
        %             y = normcdf(x, tp.mean, tp.sigma);
        %             y(y<0) = 0;
        %         case 'weibull'
        %             y = sc_in*p2p_c.weibull(tp,x);
        %         case 'power'
        %             y = sc_in*x.^tp.power;
        %         case 'exp'
        %             y = sc_in*x.^tp.k;
        %         case 'compression'
        %             y = tp.sc_in.*(tp.power.*tanh((x.*tp.sc_in)./tp.power));
        %         case 'linear'
        %             y = tp.sc_in.*x;
        %     end
        % end

        % function [p] = weibull(params, x)
        %     % [p] = Weibull(params, x)
        %     %
        %     % The Weibull function based on this equation:
        %     %
        %     % k = (-log((1-e)/(1-g)))^(1/b)
        %     % f(x) = 1 - ((1-g) * exp(-(k*x/t).^b))
        %     %
        %     % Where g is performance expected at chance, e is performance level that
        %     % defines the threshold, b is the slope of the Weibull function, and t is
        %     % the threshold
        %     %
        %     % Inputs:
        %     %   params      A structure containing the parameters of the Weibull
        %     %               function:
        %     %       b       Slope
        %     %       t       Stimulus intensity threshold as defined by 'params.e'.
        %     %               When x = 'params.t', then y = 'params.e'
        %     %       g       Performance expected at chance, proportion
        %     %       e       Threshold performance, proportion
        %     %
        %     %   x           Intensity values of the stimuli
        %     %
        %     % Output:
        %     %   p           Output of the Weibull function as a function of the
        %     %               intensity values, x
        % 
        %     % Written by G.M. Boynton - 11/13/2007
        %     % Edited by Kelly Chang - February 13, 2017
        %     % Edited by Ione Fine - February 22, 2017
        % 
        %     if ~isfield(params, 'g')
        %         params.g = 0.5;
        %     end
        %     if ~isfield(params, 'e')
        %         params.e = (0.5)^(1/3);
        %     end
        %     k = (-log((1-params.e)/(1-params.g)))^(1/params.b);
        %     p = 1 - ((1-params.g) * exp(-(k*x/params.t).^params.b));
        % end

        %% hypothetical array stuff
        % function a = Array_Sim_Location(c, a)
        %     % simulates the left visual field (vx should be negative) and right cortex (cx should be positive)
        %     % takes in parameters
        %     % c - cortical surface
        %     % a - description of the array
        %     % p  - plotting parameters
        %     %
        %     % creates various arrays:
        %     % regular_cortex  - electrodes even on visual cortex
        %     % regular_visual field - evenly spaced phosphenes in visual apce
        %     % optimal - electrode spacing based on phosphene size
        % 
        %     if nargin <2
        %         a.arrayStyle = 'optimal';
        %     end
        %     if ~isfield(a, 'eccLim'); a.eccLim = [0 32]; end % range covered by the array
        %     if ~isfield(c, 'slope') ;  c.slope = .08; end % using Keliris electrophysiology from supplementary table 1
        %     if ~isfield(c,'intercept') ;  c.intercept = .16;    end
        %     if ~isfield(a, 'rf_sz'); a.rf_sz = (c.intercept+a.eccLim*c.slope); end %  rf sizes as a funtion of eccentricity
        % 
        %     if strcmp(a.arrayStyle, 'optimal')
        %         % Packs phosphenes in visual space that vary in size as a parametric  function of eccentricity.
        %         %  Phosphene centers are then projected into
        %         % Schwartz cortical space to show how spacing of electrodes are less
        %         % densely packed near the fovea.
        % 
        %         if ~isfield(a, 'spaceFac'); a.spaceFac = 1; end % separation of the electrodes, used for the optimal array
        % 
        %         % Pack phosphenes by generating optimal spacing along the horizontal meridian by adding up
        %         % phosphene sizes
        %         xi = 0;        si = 0;      i = 1;
        %         while(xi(end)+si(end)<a.eccLim(2))
        %             si(i) = a.spaceFac*(xi(i)*c.slope+c.intercept);
        %             xi(i+1) =(xi(i)+a.spaceFac*(si(i)+c.intercept)/2)/(1-c.slope/2);
        %             i=i+1;
        %         end
        %         si(i) = a.spaceFac*(xi(i)*c.slope+c.intercept);
        % 
        %         % make concentric rings of phosphene at each spacing distance
        %         % a.vx a.vy are the positions of the phosphenes of the array in visual space
        %         ni = floor(2*pi*xi./si)';
        %         vx = 0;   vy = 0;
        %         for i=1:1:length(xi)
        %             ang = linspace(0,2*pi,ni(i)+1)';
        %             vx = [vx;xi(i).*cos(ang)];
        %             vy = [vy;xi(i).*sin(ang)];
        %         end
        %         vx = vx-.0001;   % Hack to get the foveal phosphene inside the plotting range
        %         a.nelect  = length(vx); % calculate how many electrodes we get with this spacing
        % 
        %     elseif  strcmp(a.arrayStyle, 'regular_visualfield')
        %         if ~isfield(a, 'nelect'); a.nelect =1880; end % corresponds to a spacing of 1 with fov 32 using optimal array
        %         spacing = ceil(sqrt(a.nelect));
        %         if mod(spacing,2)==1
        %             spacing = spacing +1;
        %         end
        %         tmp = linspace(-a.eccLim(2), a.eccLim(2), spacing);
        %         [vx, vy] = meshgrid(tmp, tmp); % match resolution of optimal array
        %         vx = vx(:); vy = vy(:);
        % 
        %     elseif  strcmp(a.arrayStyle, 'regular_cortex')
        %         if ~isfield(a, 'nelect'); a.nelect = 1880; end % corresponds to a spacing of 1 with fov 32 using optimal array
        % 
        %         [c_Length,~] = p2p_c.v2c_real(c, -a.eccLim(2),0);            % find axis limits
        %         [~,c_Height] = p2p_c.v2c_real(c,0,-a.eccLim(2));           % find axis limits
        %         c_Length = c_Length*1.2;     c_Height = c_Height*1.2;
        %         c_spacing = sqrt((c_Length*c_Height)/(a.nelect/2));
        % 
        %         c_Length = 0:c_spacing:c_Length;
        %         c_Height = 0:c_spacing:c_Height;
        %         c_Height = [-fliplr(c_Height(2:end)),c_Height];
        %         [cx, cy] = meshgrid(c_Length ,c_Height ); % create a grid on cortex
        % 
        %         cx = cx(:); cy = cy(:);
        %         ok = p2p_c.isValidCortex(c,cx,cy);
        %         [vx, vy] =  p2p_c.c2v_real(c,cx(ok), cy(ok));
        %         vx = cat(1, vx, -vx);vy = cat(1, vy, vy); % flip to the other side of the visual field
        %     else
        %         error('array.arrayStyle not recognized')
        %     end
        %     a.vx = vx(:); a.vy =vy(:);
        %     disp(a.arrayStyle)
        %     disp(['# electrodes = ',num2str(length(a.vx))]);
        %     a.actual_nelect = length(a.vx);
        %     a.res = c.pixpermm;
        %     a.Tfull = table(vx, vy);
        %     disp(['computed number of electrodes  = ' , num2str(length(vx))]);
        % end

        % function sz  = ecc2sz(a, ecc)
        %     % Support function for Array_Sim_Location
        %     % interpolate data to get phosphene size at locations x
        %     ecc = min(ecc, max(a.eccLim));
        %     sz = interp1(a.eccLim,a.rf_sz, ecc);
        % end

        % function Array_Sim_Plot(a, c, varargin)
        %     % Plots the array in visual space and cortical co-ordinates
        %     if nargin<3
        %         p.electrodeSize = 1; % p represents plotting parameters
        %     else
        %         p = varargin{1};
        %     end
        % 
        %     if ~isfield(p, 'electrodeSize');       p.electrodeSize =1; end
        %     if ~isfield(p, 'eccLim'); p.eccLim =[ 0, 32]; end %  [0,2.^ceil(log2(a.eccLim(2)))]; round to the nearest power of two
        %     if ~isfield(p, 'eccList'); p.eccList = [0,2.^(1:log2(p.eccLim(2)))]; end % eccentricity ring locations
        %     if ~isfield(p, 'angList'); p.angList = [90.1,135,180,225,270]*pi/180; end
        %     if ~isfield(p, 'ncol'); p.ncol = 256; end % number of colors in the colormap
        %     if ~isfield(p, 'cmap_scale'); p.cmap_scale = .25;end % how much of the colormap to use
        %     if ~isfield(p, 'symShape'); p.symShape = 'square'; end % shape of the symbols on the cortical map, 'square' or 'circle'
        %     if ~isfield(p, 'plotNums'); p.plotNums = 1:4; end
        % 
        %     % size vs eccentricity
        %     figure(p.plotNums(1));     set(gcf, 'Name' , "Phos Size vs. Ecc")
        %     ecc = linspace(p.eccLim(1),p.eccLim(2),p.ncol+1);
        %     sz = p2p_c.ecc2sz(a, ecc); % interpolate to find phosphene sizes at eccentricities defined by x
        %     idx = find(ecc<=a.eccLim(2));
        %     plot(ecc(idx), sz(idx),'-','Color',[.5,.5,.5],"LineWidth",2);
        %     set(gca,'YLim',[0,max(sz)*1.1]);    grid
        %     xlabel('Eccentricity (deg)');     ylabel('Phosphene size (deg)');
        % 
        %     % draw the phosphenes in visual coordinates, left visual field
        %     figure(p.plotNums(2)); clf; set(gcf, 'Name', 'left visual field electrodes'); hold on
        %     id = a.vx<=0;
        %     rad = sqrt(a.vx.^2+a.vy.^2);
        %     a.sz = c.slope*rad + c.intercept; % receptive field size of each electrode
        %     cmap = hsv(p.ncol);            % define colors for each electrode/phosphene
        %     col = 0.9*ones(length(id),3);  % gray
        %     scFac = p.cmap_scale*(p.ncol)/a.rf_sz(2); %scaling to get the right range in the colormap
        %     col(id,:) = cmap(min(ceil(a.sz(id)*scFac),p.ncol),:);
        % 
        %     ecc = repmat(p.eccList,p.ncol+1,1);
        %     ang = repmat(linspace(pi/2,3*pi/2,p.ncol+1)',1,length(p.eccList));
        %     gridx = ecc.*cos(ang);       gridy = ecc.*sin(ang);
        %     plot(gridx,gridy,'k-','Color',[.5,.5,.5]);
        %     ang = repmat(p.angList,p.ncol+1,1);
        %     ecc = repmat(exp(linspace(-20,log(max(p.eccList)),p.ncol+1))',1,length(p.angList));
        %     gridx = ecc.*cos(ang);    gridy = ecc.*sin(ang);
        %     cx = cos(linspace(-pi,pi,61));        cy = sin(linspace(-pi,pi,61));   % unit circle
        %     for i=1:length(a.vx)
        %         if sqrt(a.vx(i).^2+a.vy(i).^2) <a.eccLim(2)*1.4
        %             patch(a.vx(i)+a.sz(i)/2*cx,a.vy(i)+a.sz(i)/2*cy,col(i,:));
        %         end
        %     end
        %     plot(gridx,gridy,'k-','Color',[.5,.5,.5]); axis equal;
        % 
        % 
        %     m=abs(p.eccLim(2))*1.1;  set(gca,'xLim',[-m,m]);   set(gca,'YLim',[-m,m]);   axis equal
        % 
        %     % Draw the corresponding 'electrodes' in cortical space, right
        %     % hemisphere
        %     figure(p.plotNums(3));  clf;   set (gcf, 'Name', "right hemi electrodes"); hold on
        %     ecc = repmat(p.eccList,p.ncol+1,1);
        %     ang = repmat(linspace(90.1,270,p.ncol+1)'*pi/180,1,length(p.eccList));
        %     x = ecc.*cos(ang);   y = ecc.*sin(ang);
        %     [gridx,gridy] = p2p_c.v2c_real(c,x,y);
        %     plot(gridx,gridy,'k-','Color',[.5,.5,.5]);
        %     ang = repmat(p.angList,p.ncol+1,1);
        %     ecc = repmat(exp(linspace(-20,log(max(p.eccList)),p.ncol+1))',1,length(p.angList));
        %     x = ecc.*cos(ang);  y = ecc.*sin(ang);
        %     [gridx,gridy] = p2p_c.v2c_real(c,x,y);
        % 
        %     id = a.vx<0;
        %     p.z = a.vx+sqrt(-1)*a.vy;
        %     [a.cx,a.cy] = p2p_c.v2c_real(c,a.vx,a.vy);
        %     if strcmp(p.symShape, 'square')
        %         shape_x = cos(-pi/4:pi/2:5*pi/4);   shape_y = sin(-pi/4:pi/2:5*pi/4);
        %     elseif strcmp(p.symShape, 'circle')
        %         shape_x = cos(linspace(-pi,pi,61));  shape_y = sin(linspace(-pi,pi,61));
        %     end
        %     for i=1:length(a.vx)
        %         if id(i);    patch(a.cx(i)+p.electrodeSize/2*shape_x, a.cy(i)+p.electrodeSize/2*shape_y,col(i,:));  end
        %     end
        %     plot(gridx,gridy,'k-','Color',[.5,.5,.5]);
        %     [limx,limy] = p2p_c.v2c_real(c,[-.01,-.01],max(p.eccList)*1.2*[1,-1]);            % find axis limits
        %     set(gca,'XLim',[-3,limx(1)]);       set(gca,'YLim',limy*1.2);            axis equal
        % 
        %     figure(p.plotNums(4));   set(gcf, 'Name', 'Colorbar');   clf
        %     maxSz = ceil(max(a.sz)); %deg
        %     img = repmat((1:(ceil(maxSz*scFac)))',1,2);
        %     image(1,linspace(0,maxSz,size(col,1)),img);
        %     colormap(cmap)
        %     set(gca,'YDir','normal');     set(gca,'XTick',[]);      ylabel('Phosphene size (deg)'); set(gca,'FontSize',18);
        %     axis equal;       axis tight
        %     set(gcf,'PaperPosition',[1,1,1,1]);
        % end

        % function rfmaps = Array_Sim_GenerateRFmaps(a, v, c, params, range)
        %     % Creates rfmaps for optimally spaced electrodes.
        %     % takes in:
        %     % a: (optional)  -  minimally has to define the arrayStyle and
        %     % the number of electrodes in total
        %     % v:  visual field parameters (optional)
        %     % c: cortical parameters (optional)
        %     % params, containing the folowing
        %     % range: (optional) -  which values in T will be
        %     % calculated, used for manual parallelization
        %     % flag_flip: (optional) - flip RFs to save time, can be  no
        %     % flip ([ 0 0 ], up down [1 0] (top to bottom), left right [0 1] (left to right), or up down
        %     % left right from upper left visual field quadrant [ 1 1]
        %     % overwrite = 0; % skip already existing file (0) or overwrite
        %     % them (1)
        %     %
        %     % Once you've run this you can use these rfs to make movies using
        %     % Array_Sim_Movies.m
        % 
        %     v = p2p_c.define_visualmap(v);
        %     c = p2p_c.define_cortex(c);
        %     [c, v] = p2p_c.generate_corticalmap(c, v);
        %     range = range(range<=height(a.T));
        % 
        %     rfmaps = NaN(size(v.X, 1), size(v.X, 2), length(range));
        %     for i = 1: length(range)
        %         e = range(i);
        %         disp(['electrode  ', num2str(e), ' within range ', num2str(min(range)) , '-', num2str(max(range))]);
        %         v.e.x = a.T.vx(e);
        %         v.e.y = a.T.vy(e);
        %         c = p2p_c.define_electrodes(c, v);
        %         c = p2p_c.generate_ef(c);
        %         v = p2p_c.generate_corticalelectricalresponse(c, v); % create receptive field map for each electrode
        %         rfmap = mean(v.e.rfmap, 3);
        %         if ~isempty(find(isnan(rfmap(:)), 1))
        %             disp('NaNs,electrode likely outside the cortical sheet')
        %         else
        %             if params.plot
        %                 p2p_c.plotretgrid(500*(rfmap+min(rfmap(:))), v, gray(256), 2);
        %                 figure(1); p2p_c.plotcortgrid(c.e.ef, c);
        %                 drawnow;
        %                 figure(2); p2p_c.plotretgrid(255*rfmap./max(rfmap(:)), v);
        %                 drawnow;
        %             end
        % 
        %             tmp = (rfmap+params.scale(1))*params.scale(2);
        %             if max(tmp(:))>255
        %                 error('rfmap values too high for scaling to uint8 using current parameters');
        %             elseif min(tmp(:))<0
        %                 error('rfmap values too low for scaling to uint8 using current parameters');
        %             end
        %             tmp = uint8(tmp);
        %             rfmaps(:, :, i) = tmp;
        %         end
        %     end
        % end

        % function frame = Array_Sim_Movie(a, params, m)
        %     % Simulates what a movie would look like using a given array
        % 
        %     %% open input video file
        %     vid_in = VideoReader(m.filename_in);
        %     vid_dur  = vid_in.NumFrames;
        % 
        %     %% open output video file
        %     vid_out = VideoWriter(m.filename_out);
        %     vid_out.FrameRate = vid_in.FrameRate;
        %     open(vid_out);
        %     sz = size(a.rfmaps);
        % 
        %     for k = 1:vid_dur% for each frame
        %         frame = readFrame(vid_in);
        %         if k>=m.keepframes(1) && k<=m.keepframes(2)
        %             disp(['simulating frame ', num2str(k), ' out of ', num2str(vid_dur)]);
        %             if size(frame, 1)<size(frame, 2) % if the image isn't square
        %                 off = 1+ ( (size(frame, 2)-size(frame, 1))/2);
        %                 frame = frame(:, off-1:size(frame, 2)-off, :);
        %             end
        %             if ndims(frame)==3 % turn it grayscale
        %                 frame = mean(frame, 3);
        %             end
        % 
        %             in_movie=imresize(frame, [sz(1), sz(2)], "bilinear");
        %             in_movie  = in_movie./255; % scale the movie between 0 -1
        % 
        %             img = zeros(sz(1:2));
        % 
        %             for e = 1:sz(3) % for each receptive field
        %                 orig_rfmap =  (double(a.rfmaps(:, :, e))./params.scale(2))  - params.scale(1); % goes between  -1 and 1
        %                 sz_rfmap = sqrt((sum(abs(orig_rfmap(:))))./(params.pixperdeg^2));
        %                 orig_rfmap = orig_rfmap/sz_rfmap;
        %                 if isempty(find(isnan(orig_rfmap(:)), 1))
        %                     rfmap = orig_rfmap;
        %                     amp=sum(rfmap(:).*in_movie(:)); % how bright should the phosphene be, between 0 -1
        %                     img = img + (rfmap.*amp);
        %                     if params.flag_flip(1)==1
        %                         rfmap = flipud(orig_rfmap);
        %                         amp=sum(rfmap(:).*in_movie(:)); % how bright should the phosphene be, between 0 -1
        %                         img = img + (rfmap.*amp);
        %                     end
        %                     if params.flag_flip(2)==1
        %                         rfmap = fliplr(orig_rfmap);
        %                         amp=sum(rfmap(:).*in_movie(:)); % how bright should the phosphene be, between 0 -1
        %                         img = img + (rfmap.*amp);
        %                     end
        %                     if sum(params.flag_flip)==2
        %                         rfmap = fliplr(flipud(orig_rfmap));
        %                         amp=sum(rfmap(:).*in_movie(:)); % how bright should the phosphene be, between 0 -1
        %                         img = img + (rfmap.*amp);
        %                     end
        %                 end
        %             end
        %             % weirdly, because we're not really simulating the temporal part of the
        %             % model that relates electrical stimulation to brightness we don't include the nonlinearity
        %             % we assume we want the brightness to be related to the image intensity
        %             img =  img +  params.offset;
        %             figure(2)
        %             subplot(1,2,1)
        %             imagesc(in_movie(m.crop:end-m.crop, m.crop:end-m.crop)); colormap(gray(256));
        %             axis equal;    axis tight;    axis off
        %             subplot(1,2,2)
        %             image(img(m.crop:end-m.crop, m.crop:end-m.crop)); colormap(gray(256));
        %             axis equal; axis off;  axis tight
        %             drawnow
        %             frame = getframe(gca);
        %             writeVideo(vid_out,frame);
        %         end
        %     end
        %     close(vid_out);
        % end

        %% transforms
        % transforms

        % function c = define_cortical_size(c, v)
        %     vv = max(v.visfieldHeight);
        %     [ c.corticalHeight(2), ~] = p2p_c.v2c_real(c, 0, vv);
        %     vv = min(v.visfieldHeight);
        %     [c.corticalHeight(1), ~] = p2p_c.v2c_real(c, 0, vv);
        % 
        %     vv = max(v.visfieldWidth);
        %     [c.corticalWidth(2), ~] = p2p_c.v2c_real(c, vv, 0);
        %     vv =min(v.visfieldWidth);
        %     [c.corticalWidth(1), ~] = p2p_c.v2c_real(c, vv, 0);
        % end

        function [vx, vy, ok] = c2v_real(c, cx, cy)
            % takes in real x, y numbers on the cortical grid and returns
            % real x, y numbers in visual space
            [z, ok] = p2p_c.c2v_cplx(c,cx + sqrt(-1)*cy);
            vx = real(z);
            vy = imag(z);
        end

        function w = v2c_cplx(c, z)
            % takes in imaginary numbers in visual space and returns imaginary
            % positions in the cortical grid

            lvf = real(z)<0; % find points in the left visual field
            z(lvf) = z(lvf)*exp(-sqrt(-1)*pi);  % rotate lvf points 180 degrees
            w = (c.k*log(z + c.a))-c.shift; % Use the Schwartz!
            w(~lvf) = w(~lvf)*exp(-sqrt(-1)*pi);  % rotate rvf (~lvf) points back 180 degrees

            % squish
            w = real(w)+c.squish*sqrt(-1)*imag(w);
        end

        function [cx, cy] = v2c_real(c, vx, vy)
            % takes in real x, y numbers in visual space and returns real
            % x, y positions on the cortical grid
            z = p2p_c.v2c_cplx(c,vx + sqrt(-1)*vy);
            cx = real(z);
            cy = imag(z);
        end

        function [w,ok] = c2v_cplx(c, z)
            % takes in imaginary numbers in cortical space and returns
            % imaginary positions in visual space. Not all points in
            % cortical space are valid, so invalid cortical points are
            % returned as NaNs in visual space.  'ok' is a logical where 0
            % is invalid and 1 is valid.

            % unsquish
            z = real(z)+ sqrt(-1)*imag(z)/c.squish;

            lvf = real(z)>0;  % find points in left visual field (right cortex)
            z(~lvf) = z(~lvf)*exp(-sqrt(-1)*pi); % rotate rvf points 180 degres
            w = exp((z+c.shift)/c.k)-c.a;  % Undo the Schwartz!
            w(lvf) = w(lvf)*exp(-sqrt(-1)*pi);  % rotate lvf points back 180 degrees

            % set the invalid points to NaNs
            ok = p2p_c.isValidCortex(c,real(z),imag(z));
            %ok = ones(size(w));  % gmb this needs to be removed
            w(~ok) = NaN;
        end

        function ok = isValidCortex(c,cx,cy)
            % Determines which points (cx,cy) are valid in cortical space
            % Points that fall outside the C-shaped map are invalid.

            % The boundary is determined by points (0,vyb) in visual
            % cortex.  For a given y value in cortical space, cy, the
            % corresponding vyb on the boundary is:
            vyb = c.a*tan(abs(cy)/c.k);

            % The corresponding x value on the boundary in cortex, cyb, is:
            cxb = c.k*log(sqrt(vyb.^2+c.a^2))-c.shift;

            % A point in cortical space is valid if cx is greater than the
            % corresponding x point on the boundary, which is when cx>cxb,
            % and cy falls within the asyptotic range of the map (+/-
            % k*pi/2)
            ok = abs(cx)>cxb & abs(cy)<c.k*pi/2;
        end

        %% plotting functions
        % plotting functions
        % function logx2raw(base, precision)
        %     % logx2raw(base, precision)
        %     % Converts X-axis labels from log to raw values.
        %     % Inputs:
        %     %   base           Base of log transform (default: e)
        %     %   precision      Number of decimal places (default: 2)
        %     %
        %     % Example:
        %     % x = linspace(-3,0,11);
        %     % plot(log(x), log(x.^2));
        %     % logx2raw();
        %     % logy2raw(); % should be tolerant to multiple calls
        %     %
        %     % Note:
        %     % - See also: logy2raw.m
        % 
        %     % 11/17/96       gmb wrote it.
        %     % 6/6/96	     gmb added precision argument
        %     % 01/30/02       gmb updated it to use cell arrays, and to use original
        %     %                xtick values instead of converting labels. This way,
        %     %                multiple calls to this function doesn't keep converting
        %     %                the axis.
        %     % Edited by Kelly Chang - February 18, 2017
        % 
        %     %% Input Control
        % 
        %     if ~exist('base', 'var')
        %         base = exp(1);
        %     end
        % 
        %     if ~exist('precision', 'var')
        %         precision = 2;
        %     end
        % 
        %     %% Calculate Log x-axis Labels
        % 
        %     precision = sprintf('%%%2.1ff', precision*1.1);
        %     origXTick = get(gca, 'XTick'); % log x-axis labels (raw)
        %     newXTick = base.^(origXTick); % convert to raw
        %     newXLabel = arrayfun(@(x) sprintf(precision,x), newXTick, ...
        %         'UniformOutput', false); % write new x-axis labels
        %     set(gca, 'XTickLabel', newXLabel); % set x-axis labels of current graph
        % end
        % 
        % function logy2raw(base, precision)
        %     % logy2raw(base, precision)
        %     %
        %     % Converts Y-axis labels from log to raw values.
        %     %
        %     % Inputs:
        %     %   base           Base of log transform (default: e)
        %     %   precision      Number of decimal places (default: 2)
        %     %
        %     % Example:
        %     % x = linspace(-3,0,11);
        %     % plot(log(x), log(x.^2));
        %     % logx2raw();
        %     % logy2raw(); % should be tolerant to multiple calls
        %     %
        %     % Note:
        %     % - See also: logx2raw.m
        % 
        %     % 11/17/96       gmb wrote it.
        %     % 6/6/96	     gmb added precision argument
        %     % 01/30/02       gmb updated it to use cell arrays, and to use original
        %     %                xtick values instead of converting labels. This way,
        %     %                multiple calls to this function doesn't keep converting
        %     %                the axis.
        %     % Edited by Kelly Chang - February 18, 2017
        % 
        %     %% Input Control
        % 
        %     if ~exist('base', 'var')
        %         base = exp(1);
        %     end
        % 
        %     if ~exist('precision', 'var')
        %         precision = 2;
        %     end
        % 
        %     %% Calculate Log x-axis Labels
        % 
        %     precision = sprintf('%%%2.1ff', precision*1.1);
        %     origYTick = get(gca, 'YTick'); % log y-axis labels (raw)
        %     newYTick = base.^(origYTick); % convert to raw
        %     newYLabel = arrayfun(@(x) sprintf(precision,x), newYTick, ...
        %         'UniformOutput', false); % write new y-axis labels
        %     set(gca, 'YTickLabel', newYLabel); % set y-axis labels of current graph
        % end

        function plotcortgrid(img, c,  varargin)
            % plotcortgrid(img, c, cmap,figNum, evalstr)
            % takes as input:
            %   cortical image you want to plot, such as c.e.ef
            %   the structure c that defines the cortical surface
            % optional arguments:
            %   colormap, figure number and a string to evaluate
            %  (e.g. 'title("electric field")' or 'subplot(1, 2,1)';
            %
            % commented 12/7/2024 IF

            if nargin<3 || isempty(varargin{1});  cmap = gray(256);   else; cmap = varargin{1}; end
            if nargin<4 || isempty(varargin{2});  figNum = 1;  else; figNum = varargin{2}; end
            if nargin<5 || isempty(varargin{3});  evalstr = '';      else;  evalstr = [varargin{3}]; end

            if isfield(c,'cropPix')
                img(isnan(c.cropPix)) = NaN;
                img = img + 2;
                cmap = [0,0,0;cmap];
            end

            eval(evalstr); colormap(cmap);
            if ~isempty(img)
                if min(img<2)
                    image(c.x, c.y, img+20); hold on
                else
                    image(c.x, c.y, img); hold on
                end
            end
            xlabel('mm'); ylabel('mm')
            set(gca,'YDir','normal');
            if mean(c.cortexLength)<=0
                plot(c.v.gridAng, '-', 'Color', c.gridColor); hold on
                plot(c.v.gridEcc, '-', 'Color', c.gridColor);
            end
            if mean(c.cortexLength)>=0 %
                % plotting both hemispheres symmetrically
                plot(-c.v.gridAng, '-', 'Color', c.gridColor);
                plot(-c.v.gridEcc, '-', 'Color', c.gridColor);
            end

            axis equal;  axis tight
            set(gca,'XLim',[min(c.x(:)),max(c.x(:))]);
            set(gca,'YLim',[min(c.y(:)),max(c.y(:))]);
            drawnow;
        end

        function plotretgrid(img, v, varargin)
            % plotretgrid(img, c)
            % plotretgrid(img, c, cmap,figNum, evalstr)
            % takes as input:
            %   retinal image
            %   the structure v that defines the retinal surface
            % optional arguments:
            %   colormap, figure number and a string to evaluate
            %  (e.g. 'title("visual field")' or 'subplot(1, 2, 1)';
            %
            % commented 12/7/2024 IF

            if nargin<3 || isempty(varargin{1});  cmap = gray(256);   else; cmap = varargin{1}; end
            if nargin<4 || isempty(varargin{2});  figNum = 1;        else; figNum = varargin{2}; end
            if nargin<5 || isempty(varargin{3});  evalstr = '';      else; evalstr = [varargin{3}]; end

            figure(figNum); hold on
            eval(evalstr);
            image(v.x, v.y, img); hold on

            colormap(cmap); set(gca,'YDir','normal');

            plot(v.zAng,'-','Color', v.gridColor); plot(v.zEcc,'-','Color', v.gridColor);
            plot(-v.zAng,'-','Color', v.gridColor); plot(-v.zEcc,'-','Color', v.gridColor);

            axis equal;  axis tight
            xlabel('degrees'); ylabel('degrees')
            set(gca,'XLim',[min(v.x(:)),max(v.x(:))]);
            set(gca,'YLim',[min(v.y(:)),max(v.y(:))]);
            drawnow;
        end

        function p = fit_ellipse_to_phosphene(img,v)
            M00 = sum(sum(img));
            M10 = sum(sum(v.X.*img));           M01 = sum(sum(v.Y.*img));         M11 = sum(sum(v.X.*v.Y.*img));
            M20 = sum(sum(v.X.^2.*img));          M02 = sum(sum(v.Y.^2.*img));
            p.x0 = M10/M00;         p.y0 = M01/M00;
            mu20 = M20/M00 - p.x0^2;      mu02 = M02/M00 - p.y0^2;             mu11 = M11/M00 - p.x0*p.y0;
            a = (mu20+mu02)/2;         b = .5*sqrt(4*mu11^2+(mu20-mu02)^2);
            lambda_1 = a+b;      lambda_2 = a-b;
            p.theta = -.5*atan2(2*mu11,mu20-mu02);
            p.sigma_x = 2*sqrt(lambda_1);        p.sigma_y = 2*sqrt(lambda_2);
        end

        % function fillSymbols(h,colList)
        %     if ~exist('h', 'var');     h = get(gca,'Children');    end
        %     for i=1:length(h)
        %         if ~exist('colList','var');       col = get(h(i),'Color');
        %         else
        %             if iscell(colList); col = colList{i};
        %             else;     col = colList(i,:);  end
        %         end
        %         set(h(i),'MarkerFaceColor',col);
        %     end
        % end

        % function draw_ellipse(trl, figNum, spstr, varargin)
        %     figure(figNum); hold on
        %     eval(spstr);
        %     if nargin >3
        %         eye = varargin{1};
        %     else
        %         eye = 1;
        %     end
        % 
        %     if nargin>4
        %         lineColor = varargin{2};
        %     else
        %         lineColor = 'g';
        %     end
        %     theta = linspace(-pi,pi,101);
        % 
        %     for e=1:length(eye)
        %         r = sqrt( (trl.ellipse(eye(e)).sigma_x*trl.ellipse(eye(e)).sigma_y)^2./ ...
        %             (trl.ellipse(eye(e)).sigma_y^2*cos(theta).^2 + trl.ellipse(eye(e)).sigma_x^2*sin(theta).^2));
        %         x = trl.ellipse(eye(e)).x+r.*cos(theta-trl.ellipse(eye(e)).theta);
        %         y = trl.ellipse(eye(e)).y+r.*sin(theta-trl.ellipse(eye(e)).theta);
        %         plot(x,y,'-','LineWidth',1,'Color',lineColor);
        %     end
        % end

        % function  err = (p,trl, v)fit_phosphene
        %     % finds the standard deviation and mean of the best fitting
        %     % 2D Gaussian, assumes circularity.
        %     phos = mean(trl.max_phosphene, 3); phos = phos./max(phos(:));% average across the two eyes
        %     pred = p2p_c.Gauss(p, v);  pred = pred./max(pred(:));
        % 
        %     err = sum((phos(:)-pred(:)).^2);
        %     if p.sigma<0; err =err +  abs(p.sigma)*10.^6; end
        % end

        % function out = Gauss(p,v)
        %     out = exp(-((v.X-p.x).^2+(v.Y-p.y).^2)/(2*p.sigma^2));
        % end

%         function [err,predcx,predcy,cx,cy] = fitElectrodeGrid(p,vx,vy)
%             % Gives a fit of a projected set of phosphenes onto a 6x4
%             % electrode grid by first projecting the phosphenes using the
%             % v2c mapping function and then moving and rotating the grid.
% 
%             p.shift = p.k*log(p.a);  % This really should be built in
% 
%             % Project phosphenes into cortex (including squishing)
%             [cx,cy] = p2p_Beauchamp.v2c_real(p, vx, vy);
% 
%             % Creat a 6x4 array
%             [x0,y0] = meshgrid(-1.5:1.5,-2.5:2.5);
%             x0 =x0(:)';
%             y0 =y0(:)';
% 
%             % Expand by dx
%             x0 = x0*p.dx;
%             y0 = y0*p.dx;
% 
%             % Rotate by ang and shift by (xc,yc)
%             rot = [cos(p.ang) sin(p.ang);-sin(p.ang) cos(p.ang)];
% 
%             M = rot*[x0;y0] + repmat([p.xc;p.yc],1,24);
%             predcx = M(1,:);
%             predcy = M(2,:);
% 
%             % Compare predicted to projected
%             err = sum((predcx-cx).^2 + (predcy-cy).^2);
%         end
    end
end

