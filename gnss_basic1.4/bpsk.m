function [timetags, baseband_signal, chipindex] = ...
    bpsk ( timetags, nsamples, prncodes, prnflags);
%
%  Generates a synthetic GPS signals at the same sample rate as the data
%
%  sample_rate = Hz
%  nsamples = total number of samples to generate
%  prncodes - 1023 X 32 array of all PRN codes, 
%  prnflags - Which of the PRN codes you want to generate a baseband
%   signal for.
%
%  timetags - sample times for the generated signal - in seconds.
%  baseband_signal - BPSK signal, modulated from -1 to +1.
%
Tc = 9.775E-7;
%timetags = [0:1:nsamples]' / sample_rate; % in sec.
chipindex = floor((1/Tc) * timetags ); % in chips
chipindex = mod(chipindex, 1022) + 1;  % allow the code to "roll-over"
baseband_signal = -2*prncodes(chipindex,prnflags) + 1;

return