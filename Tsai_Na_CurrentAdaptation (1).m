% Tsai_Na_CurrentAdaptation

% we believe silencing is based on slow recovery from adaptation in sodium
% channels
% https://iopscience.iop.org/article/10.1088/1741-2560/8/6/066007/pdf
% Frequency-dependent reduction of voltage-gated sodium current modulates retinal ganglion cell response rate to electrical stimulation
%  David Tsai et al 2011 J. Neural Eng. 8 066007


rng(1)
clc; clear
close all
addpath('C:\Users\Ione Fine\Documents\code\toolboxes\UWToolbox');
TSAI = readtable('C:\Users\Ione Fine\Documents\code\p2p-cortical\datasets\Tsai_NaCurrent_2011.csv');

% define parameters

trl.amp = 5000;
trl.pw =  .1 * 10^-3;
trl.dur = 0.13;
figure(1)
clrList = [1 0 0 ; .5 .0 .5; .3 0 1];

tp.saturation_model = 'linear';
tp.spikemodel = 'integratefire';
tp.gammaflag = 1;

tp.Na_recovery = 2.6457;
tp.Na_strength =3.1345e-04;
tp = p2p_c.define_temporalparameters(tp);
plt_index = 1:30000;
freqList =  [50 100 200];


%% load Tsai data showing sodium current decline with repeated pulses
for f = 1:3
    for np = 1:6
        ind = find(TSAI.NumberPulses==np & TSAI.Freq == freqList(f));
        Na_current(f,np) = TSAI.NaCurrent(ind);
        trl.freq  = TSAI.Freq(ind);
    end
end

freeList = {'Na_recovery', 'Na_strength'};
tp = optUW.fit('fit_Tsai_Na_current', tp, freeList,  trl, freqList, Na_current);

pred_Na_current = calc_Tsai_Na_current(tp, trl, freqList);

for f = 1:3
    trl.ip = 0;
    trl = p2p_c.define_trial (tp, trl);
    trl = p2p_c.spike_model(tp, trl, 20);
    figure(f); clf
    subplot(2, 1, 1)
    % plot the pulse train
    plot(trl.t(plt_index),trl.pt(plt_index)./max(trl.pt), 'k'); hold on
    plot(trl.t(plt_index),trl.Na_adapt(plt_index), 'b'); hold on
    set(gca, 'YLim' , [-1.1 1.1] )
    % plot predicted and actual Na current over time
    subplot(2,1,2)
    plot(1:6, Na_current(f,:), 'o', "MarkerFaceColor", clrList(f,:), "MarkerEdgeColor", clrList(f,:)); hold on
    plot(1:6, pred_Na_current(f,:), '-', "Color", clrList(f,:)); hold on
    set(gca, 'YLim', [0 1])
end
%     cutoff = 6*(1/trl.freq);
% 
%     subplot(3,1,2)
%     plot(trl.t(plt_index), trl.Na_adapt(plt_index)); hold on
% 
% 
% 
% for f = 1:3
%     trl.freq = freqList(f);
%     trl = p2p_c.define_trial (tp, trl);
%     trl.pt(trl.t>cutoff) = 0; % blank out some pulses
%     trl = p2p_c.spike_model(tp, trl, 20);
% 
%     for np = 1:6
%         mask = zeros(length(trl.pt), 1);
%         mask(trl.t>(((np-1)/trl.freq)-.0025) & trl.t<=((np/trl.freq)+.0025)) = 1;
%         adapt(f, np) = max(trl.Na_adapt.*mask);
%     end
% end
% legend({'50', '100', '200'})
% 
% 
% return
% for t = 1:size(TSAI, 1)
%     trl.freq  = TSAI.Freq(t);
%     trl = p2p_c.define_trial (tp, trl);
%     figure(2); clf
%     subplot(3, 1, 2)
%     plot(trl.t(plt_index), sum(trl.indexFire(plt_index, :), 2), 'b');                                         
%     disp([' number of spikes  = ', num2str(sum(trl.indexFire(:)))]);
%     subplot(3, 1, 3)
%     plot(trl.t(plt_index),trl.resp(plt_index), 'b'); hold on
%     pause
% end