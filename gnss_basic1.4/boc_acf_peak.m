function [tau, peak_ts, peak_tc, RY_peak] = boc_acf(fs, fc)
%
% Infinite BW ACF for BOC(alpha, beta) modulated signal.
%  peak_ts = locations in spreading period ts of the peaks
%  peak_tc = locations in the code chipping period tc of the peaks
%  RY_peak = ACF peak values at the points in peak_ts or peak_tc
n = 2 * fs/fc
ell = -n:1:n;
peak_ts = ell';
peak_tc = peak_ts/n;
RY_peak = ((-1).^ell ) .* ( n - abs(ell))/n;

peak_tc = [-10; peak_tc; 10];
RY_peak = [0; RY_peak'; 0];

fc0 = 1.023E6;
tc = 1/(fc * fc0);
tau = tc * peak_tc;

return