function [R2, R] = boc_acf( tau, fs, fc)
% 
% ACF for BOC(fs, fc) modulated signals
%  tau is in units of TC   (tau =  tau_sec / Tc_sec) to stay consistent
%  with the other ACF functions. 
% fs, fc are the multiple of 1.023 MHz
%
%  See Betz, 2001 Navigation paper.

n = 2 * fs/fc
ell = -n:1:n;   % index of peaks (2n-1 total)
peak_ts = ell';  % Location of peaks, in units of Ts
peak_tc = peak_ts/n;  % location of peaks in units of Tc.
RY_peak = ((-1).^ell ) .* ( n - abs(ell))/n; % value of peaks

% Values of the ACF at each peak.
peak_tc = [min(tau); min(peak_tc) - 1/n; peak_tc; max(peak_tc) + 1/n; max(tau)];
RY_peak = [0; 0; RY_peak'; 0; 0];

% Interpolated to the lags, tau

R = interp1(peak_tc, RY_peak, tau);
R2 = R.^2;

return