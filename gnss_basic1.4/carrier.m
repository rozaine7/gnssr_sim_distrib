function [I,Q] = ...
    carrier( sample_rate, baseband_signal, fcarrier)
%
% Generates a simulated down-converted GPS signal by modulating the
% sampled PN codes (baseband_signal) onto the inphase(I) and
% quadrature(Q) carrier at a frequency of fcarrier
%
% fcarrier and sample_rate are in Hz.
%
[nsamples, ncodes] = size(baseband_signal);
phase = ([0:nsamples-1]' * (fcarrier/sample_rate)); % in cycles
I = baseband_signal .* cos( 2 * pi * phase );
Q = baseband_signal .* sin( 2 * pi * phase );

return

