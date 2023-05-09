function [f, Sf] = boc_acf(f1, f2, num, den, frange)
%
% Infinite BW Spectrum for MBOC(fs, fc, num/den) modulated signal.
%  frange = range of frequencies, in MHz
%  fs = Frequency, units of 1.023MHz
%  fc = Frequency, units of 1.023MHz
%  Sf - power spectrum
% The nomenclature is:  BOC(f2,f2) and BOC(f1, f2) are linearly combined.

f = linspace(-frange, frange, 1000);
fc0 = 1.023e6;  % Hz - just so the units turn out OK.
fc1MHz = f2;
fs1MHz = f2;
fc2MHz = f2;
fs2MHz = f1;

[fdummy, Sf1] = boc_spec(fs1MHz, fc1MHz, frange);
[fdummy, Sf2] = boc_spec(fs2MHz, fc2MHz, frange);

Sf = (den-num)/den * Sf1 + (num/den) * Sf2;

return