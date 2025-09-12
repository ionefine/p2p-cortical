function pred_Na_current = calc_Tsai_Na_current(tp, trl,freqList)

for f = 1:3
    trl.freq = freqList(f);
    trl = p2p_c.define_trial (tp, trl);
    trl = p2p_c.spike_model(tp, trl, 20);

    ss = .5/freqList(f);

    figure(f); clf
    for np = 1:6
        mask = trl.pt;
        end_cutoff = trl.ip+ss+(ss*2*(np-1));
        st_cutoff = ss+(ss*2*(np-2));
        mask(trl.t>end_cutoff) = 0; % blank out some pulses
        mask(trl.t>end_cutoff) = 0; % blank out some pulses
        mask(trl.t<st_cutoff) = 0; % blank out some pulses
        mask = mask./max(mask(:));
        mask(trl.pt<0) = 0;
        pred_Na_current(f, np) = max(trl.Na_adapt.*trl.pt);
        % subplot(6,1, np)
        % plot(mask(1:120000)/max(mask(:))); hold on
        % plot(trl.Na_adapt(1:120000))
        pred_Na_current(f, np) = mean(trl.Na_adapt(find(mask)));
    end
end
