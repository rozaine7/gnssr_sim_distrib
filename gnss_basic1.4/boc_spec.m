function [f, Sf] = boc_acf(fs, fc, frange)
%
% Infinite BW Spectrum for BOC(fs, fc) modulated signal.
%  frange = range of frequencies, in MHz
%  fs = Frequency, units of 1.023MHz
%  fc = Frequency, units of 1.023MHz
%  Sf - power spectrum
% fixed some unit problems 4/15/21 - JLG - everything calculated in Hz and
% s  Sf should be in W/Hz
%
n = 2 * fs/fc

f = linspace(-frange*1e6, frange*1e6, 1000); % HZ
fc0 = 1.023e6;  % Hz - just so the units turn out OK.
fcHz = fc * fc0;
fsHz = fs * fc0;
Ts = 1./(2*fsHz); % Ts is subcarrier half-period
%fprime = f / fsMHz;

if (floor(n/2) == n/2)   % n even
  Sf = (1/(n*Ts)) .* ( sin( pi .* f * Ts) .* sin ( n .* pi .* f .*Ts)./(pi .* f .* cos(pi .* f.* Ts))).^2;
else  % n odd
 Sf = (1/(n*Ts)) .* ( sin( pi .* f * Ts) .* cos ( n .* pi .* f .*Ts)./(pi .* f .* cos(pi .* f.* Ts))).^2;
end

return