function [steady_state_dyn, sigma_noise] = dll_err ( ...
          BL, Ti, d, SN0, los_dyn, flag)
%
%  Thermal (random) noise error + dynamics error for an E-L DLL
%
%  Flag used to set loop order 
%
%  BL = close loop BW in Hz
%  Ti = Predetection Integration time (sec).
%  SN0 = Carrier to noise ratio (non-dim, linear scale)
%  los_dyn = line of sight dynmaics, for 2nd order DLL = Acceleration
%  (dr^2/d^2t) in m/sec^2
%
if (flag == 2)
   omega0 = 1.89 * BL; % for critically damped loop
   los_dyn_chips = abs(los_dyn) / 293; % dynamics in chips/sec^2

   steady_state_dyn = (1./omega0.^2) * los_dyn_chips;  %chips
elseif (flag == 3)
   omega0 = 1.2 * BL; % for critically damped loop
   los_dyn_chips = abs(los_dyn) / 293; % dynamics in chips/sec^2

   steady_state_dyn = (1/omega0.^3) * los_dyn_chips;  %chips
else
   steady_state_dyn = 0;
end

sigma2_noise = (BL * d)./(2 * SN0) .* (1 + 2./((2-d)*SN0*Ti));

sigma_noise = sqrt(sigma2_noise);

return

