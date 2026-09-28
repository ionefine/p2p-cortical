% Simulate_Elche_Drawings_Main
% Written by IF 3/7/2025
% Commented by ES 7/25/25
% functions added by ES September 2025
%
% takes in a patient drawing and finds percepts that looks similar
% then percepts are combines to better match images
% first combinations with best correlation are found
% then random combinations are found
%
% note if you want to run this whole thing you might want to go through and
% changed the names of thing you are saving and loading as to not overwrite
% previous things/cause a lot of things are saved and that way you know
%
% also if a function works until it hits a certain drawing try clearing
% everything unneccessary and running it just for that drawing


clear
close all

% sets up the drawing and helper func folders
vbl = ed.setup(computer); % or 'Z' computer
vbl.cleanstartflag = 0;


% save cropped version of both original and edited versilocon
% has to be done once at least
% ed.crop_drawings(vbl);
p_draw = ed.load_patient_drawings(vbl);

%% define cortical model
[c,v, trl, tp] = ed.define_cortical_model(vbl);
if 1

%% generate + save "best" phosphenes
p_draw = ed.simulate_drawings(c, v, trl, tp, p_draw, vbl);
ed.save_simulated_drawings(p_draw, vbl);

%ed.visualize_singletons(vbl,7);

%% combine the drawings
ed.combine_sim_draws(vbl);
for d = 1:length(vbl.dirList)
    ed.visualize_combinations(vbl, 7);
end
end
%% create the random versions

c = rmfield(c, 'e');v = rmfield(v, 'e');
c = rmfield(c, 'x'); c = rmfield(c, 'y');v = rmfield(v, 'x'); v = rmfield(v, 'y');
c = rmfield(c, 'X'); c = rmfield(c, 'Y');
c = rmfield(c, 'v'); c = rmfield(c, 'cropPix');

parpool(15)
for i=25:100    
    parfor rep = (i-1)*15+1:(i-1)*15+15
        if vbl.debugflag
            filename = fullfile(vbl.datadir, vbl.dirList(end).name, 'random_models', [vbl.dirList(end).name, vbl.fileidstr, '_', num2str(rep), '_debug.mat']);
        else
            filename= fullfile(vbl.datadir, vbl.dirList(1).name, 'random_models', [vbl.dirList(end).name, vbl.fileidstr, '_', num2str(rep), '.mat']);
        end
        if exist(filename, 'file')
            disp([num2str(rep), ' exists']);
        else
            ed.simulate_save_rand_drawings(rep, c, v, trl, tp,  vbl);
        end
    end
    delete(gcp);
end
disp('done with generating random images')


%% create random combinations
ed.combine_random_models(vbl);


%% analysis rand vs. realend
% these functions work for one drawing at a time
% just rerun this part to look at different drawings
for d = 1:length(vbl.dirList)
    %ed.show_best_combos(vbl, d); % show top real and rand combos images
    %ed.plot_corr_histograms(vbl, d)
    ed.compare_real_vs_rand_stats(vbl, d)
end
ed.plot_corr_histograms(vbl, d) % plot best real corr vs best rand per rep corr on a histogram
ed.compare_real_vs_rand_stats(vbl, d) % look at the percentile where the top real lands in the distribtion of top rands per rep


