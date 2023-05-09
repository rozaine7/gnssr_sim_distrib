function [xquant, discrval] = quantize( x, nbits, minval, maxval )
%
% Crude implementation of quantizing a floating point signal
%
nlevels = 2^nbits;

step = (maxval-minval)/nlevels;

discrval = floor((x-minval)/step);

xquant = minval + discrval*step;


end

