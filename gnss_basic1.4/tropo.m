function [TD, ZTW, ZTD, MD, MW] = tropo( PR, TC, RH, H, phi, EL, ...
              modelflag, mapflag)
%
% Tropospheric models:
%  modelflag = 'h' -> Hopfield
%              's' -> Saastamoinen
%
%  mapflag   = 1 -> 1/sinE
%            = 2 -> 1/sqrt(1-(cosE/1.001)^2);
%            = 3 -> MD and MW separate
% Inputs = PR -> Pressure (mbar) [matrix]
%          TC -> Temp in Celsius [matrix]
%          RH -> Relative humidity. [matrix]
%          H -> Orthometric ht (m) [scalar]
%          phi -> latitude       [scalar]
%
TK = TC + 273.15;
hkm = H / 1000;
phirad = phi * pi/180;
ELrad = EL * pi/180;
%
% Zenith delay models
%
if ( modelflag == 's')
   fprintf(' Saastamoinen Model \n');
   e0 = 6.108 * RH .* exp( (17.15*TK -4684)./(TK - 38.45));
   ZTD = 0.0022777*( 1 + 0.0026 * cos( 2*phirad) + 0.00028 * hkm) * PR;
   ZTW = 0.0022777*(1255./TK + 0.05).*e0;
elseif (modelflag == 'h')
   fprintf(' Hopfield Model \n');
   hd = 40136 + 148.72 * TC; %m - from H-W
   hw = 12e3; %m
   e0 = 6.108 * RH .* exp( (17.15*TK -4684)./(TK - 38.45));
   ZTD = (77.64E-6 * PR .* hd) ./ (TK * 5);
   ZTW = (0.373 * e0 * hw) ./ (TK.^2 * 5);
else
   fprintf(' Invalid model type \n')
   return
end

if (mapflag == 1)
   MW = 1./sin(ELrad);
   MD = 1./sin(ELrad);
elseif (mapflag == 2)
   MW = 1./sqrt(1 - (cos(ELrad)/1.001).^2);
   MD = 1./sqrt(1 - (cos(ELrad)/1.001).^2);
elseif (mapflag == 3)
   MW = 1./( sin(ELrad) + 0.00035./(tan(ELrad) + 0.017));
   MD = 1./( sin(ELrad) + 0.00143./(tan(ELrad) + 0.0445));
else
   fprintf('Invalid mapping function \n');
   return
end

TD = ZTW.*MW + ZTD.*MD;