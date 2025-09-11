function pred_Na_current = calc_Tsai_Na_current(tp, trl,freqList)

for f = 1:3
    trl.freq = freqList(f);
    trl = p2p_c.define_trial (tp, trl);
    trl = p2p_c.spike_model(tp, trl, 20);

    for np = 1:6
        mask = zeros(1,length(trl.pt));
        mask(trl.t>(((np-1)/trl.freq)+.0025) & trl.t<=((np/trl.freq)+.0025)) = 1;
        mask = mask.*max(trl.pt/max(trl.pt), 0);
        pred_Na_current(f, np) = max(trl.Na_adapt.*mask);
    end
end
