function [lagc, corr_pwr] = ...
    code_search( sample_rate, if_freq, doppler_freq, data_set, prncodes, prnflags);
%
%  Computes the correlation between locally generated signal at the 
%  down-converted carrier frequency of carrier_freq (Hz), for the
%  satellite PRNs identified in prnflags.
%
%  The data is in the column-vector data_set, sampled at the rate of 
%  sample_rate (Hz).
%
%  The cross-correlation is computed at lag intervals of 1/2 chip
%
nsamples_data = size(data_set,1);
Tc = 9.775e-7; % sec
samples_per_chip = sample_rate  * Tc;
carrier_freq = if_freq + doppler_freq;
L1 = 1575.42e6;

timetags = [0:nsamples_data-1] * (1-doppler_freq/L1)/sample_rate;
[timetags, baseband_signal, cid] = ...
    bpsk( timetags, nsamples_data-1, prncodes, prnflags);
[I,Q] = carrier( sample_rate, baseband_signal, carrier_freq);
                                                 % steps
[Icorr,lags] = ...
    xcorr(data_set, I, floor(1023 * 0.5 * samples_per_chip), 'unbiased');
Qcorr = xcorr(data_set, Q, floor(1023  * 0.5 * samples_per_chip), 'unbiased');

corr_pwr = Icorr.^2 + Qcorr.^2;

lagc = lags'/ samples_per_chip;

return
