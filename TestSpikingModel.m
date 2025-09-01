% TestSpikingModel
% compares spiking model to various data and the convolve model
% written ES & IF 2025
rng(1)
clc;
close all

% temporal parameters
tp = p2p_c.define_temporalparameters();

v.drawthr = 1;

% define parameters
freq  = [50 100 200];
freq  = [50  200];
amp = 500*ones(size(freq));
dur =  70*10^-3*ones(size(amp));
pw =  .1 * 10^-3*ones(size(amp));

% create table to hold parameter combos
Tsim = table(amp', pw', freq', dur');
Tsim.Properties.VariableNames = {'amp', 'pw','freq','dur'};

% generate trials with different amplitudes but fixed duration, freq, pw
all_trl = p2p_c.loop_model(tp, Tsim);
plt_index = 1:length(all_trl(1).pt)/3;
for f = 1
    disp(['freq - ', num2str(f)])
    for npulses = 1
        disp(['pulse - ', num2str(npulses)])
        trl = all_trl(f);
        trl.pt = trl.pt(plt_index);
         trl.t = trl.t(plt_index);
        ss = .5/freq(f);
        cutoff = ss+(ss*2*(npulses-1));
        trl.pt(trl.t>cutoff) = 0; % blank out some pulses

        tp.saturation_model = 'linear';
        tp.spikemodel = 'integratefire';
        trlS = p2p_c.spike_model(tp, trl, 20);
        respS(f, npulses) = max(trlS.resp);
        plot(trlS.t(plt_index),trlS.pt(plt_index)/50, 'k'); hold on
        plot(trlS.t(plt_index),trlS.resp(plt_index), 'b'); hold on

        tp.saturation_model = 'compression';
        tp.spikemodel = 'convolve';
        tp.gammaflag = 1;
        trlC = p2p_c.spike_model(tp, trl);
        respC(f, npulses) = max(trlC.resp);
        plot(trlC.t(plt_index),trlC.resp(plt_index), 'r')
    end
end

clist = ['r','g', 'b'];
for f = 1:length(freq)
    plot(1:6, respS(f,:), 'Color', clist(f)); hold on
    xlabel('number of pulses');
    ylabel('spiking response')
end
