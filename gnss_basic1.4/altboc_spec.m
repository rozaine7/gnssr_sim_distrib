function [f, Sf, Sfcm] = boc_spec(fs, fc, frange)
%
% Infinite BW Spectrum for altBOC(fs, fc) modulated signal.
%  frange = range of frequencies, in MHz
%  fs = Frequency, units of 1.023MHz
%  fc = Frequency, units of 1.023MHz
%  Sf - power spectrum
n = 2 * fs/fc
f = linspace(-frange, frange, 1000);
fc0 = 1.023e6;  % Hz - just so the units turn out OK.
fcMHz = fc * fc0;
fsMHz = fs * fc0;
fprime = f / fcMHz;

Ts = 1/fsMHz;
Tc = 1/fcMHz;

if (floor(n/2) == n/2)   % n even
 fprintf(' n even not implemented yet - only valid for odd n \n')
else
  fprintf(' n = %4i \n', n);
  Sf = (8 ./ (Tc * pi^2 * f.^2)).* ...
       (cos(pi * f * Tc)./cos(pi * f * Tc/n)).^2 .* ...
       ( 1 - cos( pi * f * Tc/n));
  
  Sfcm = (4 ./ ( pi^2 * f.^2 * Tc)) .* ...
         (cos( pi * Tc * f)./cos(pi * f * Tc /n)).^2 .* ...
         ( (cos( pi * f * Ts/2)).^2 - cos( pi * f * Ts /2 ) - ...
            2 * cos( pi * f * Ts/2).*cos(pi * f * Ts/4) + 2 );
end

return