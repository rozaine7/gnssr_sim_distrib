function [xquant, discrval] = quant( x, nbits, minval, maxval )
%
% Crude implementation of quantizing a floating point signal
%  nbits <=0 will return the floating point value. 
%
if nbits > 0
 nlevels = 2^nbits;

 step = (maxval-minval)/nlevels;

 discrval = floor((x-minval)/step);

 xquant = minval + discrval*step;

else
    xquant = x;
    discrval =0;

end

