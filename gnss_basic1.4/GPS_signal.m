function [I, Q, timetags, chipindex] = GPS_signal ( sample_rate, nsamples, ...
			carrier_DCO, code_DCO, prncodes, prnflags);
%
%  Generates a synthetic GPS signals at the same sample rate as the data
%
%  sample_rate = Hz
%  nsamples = total number of samples to generate
%  carrier_DCO = post-downconverison carrier (Hz)
%  code_DCO = chips/sec 
%
%  This is most general
%
timetags = [0:1:nsamples]' / sample_rate;
chipindex = mod( floor( code_DCO * timetags), 1023) + 1;
size(timetags)
size(chipindex)
I = prncodes(chipindex,prnflags) .* ...
    (cos( 2 * pi * carrier_DCO * timetags ) * ones(1,size(prnflags,2)));
%Q = prncodes(chipindex,:) .* ...
%    (sin( 2 * pi * carrier_DCO * timetags ) * ones(1,32));
return