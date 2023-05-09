function [t, x, xmod ] = boc_mod(fs, fc, bits)
%
% generate a plot of the boc-code modulation of the bit sequence in bits.
%
%  fs = subcarrier frequency in units of 1.023 MHz
%  fc = chipping frequency in units of 1.023 MHz
%  
nbits = max(size(bits));
fc0 = 1.023E6;
fsHz = fs * fc0;
fcHz = fc * fc0;

Tc = 1/fcHz;
Ts = 1/(2*fsHz);

n = 2 * fs/fc

t = 0:Ts/20:nbits*Tc-eps;

modindex = floor(t/Tc)+1;

xmod = 1- 2* bits(modindex);

subcarrier = (-1).^[0:n-1];

subcarrier_index = mod(floor(t/Ts),n) +1 ;

xsub = subcarrier(subcarrier_index);

if floor(n/2) == n/2  % n even
  x = xmod .* xsub;
else
  x = (-1).^(modindex - 1) .* xmod .* xsub;
end

  