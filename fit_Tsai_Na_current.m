function err = fit_Na_current(tp, trl, freqList, Na_current)


pred_Na_current = calc_Tsai_Na_current(tp, trl, freqList);

err = sum((pred_Na_current(:)-Na_current(:)).^2)