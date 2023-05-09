function S = elfouhaily1d(k, U10, inverse_wave_age)
%
% Generates the omnidirectional Elfouhaily height spectrum
% as a function of wavenumber, k in units of rad/m.  
%
% The input parameters are U10=wind speed at 10 m
% elevation (m/s) and inverse_wave_age (nondimensional).
%
%
% Long wave curvature spectrum. (section 5.2.1)
%
% Corrected based on errors found in UPC class (6/07).  Equations should
% now match those in the Elfouhaily 1997 paper.
%
g = 9.81; % m/sec^2
cm = 0.23; %m/s
u_star = 0.023 * U10^(1.23);
km = 370; % rad/m - from text on p 15789 - not clear if this is constant.
k0 = g / (U10^2);
cp = U10/inverse_wave_age;
c = sqrt(g * (1 + (k/km).^2) ./ k);
alpha_p = 6E-3 * sqrt(inverse_wave_age); % eq. (34)
sigma = 0.08 * (1 + 4 / (inverse_wave_age^3));
kp = k0 * inverse_wave_age^2;
if ((inverse_wave_age < 1) & (inverse_wave_age >= 0.84))
  littlegamma = 1.7;
elseif ((inverse_wave_age >= 1) & (inverse_wave_age < 5 ))
  littlegamma = 1.7 + 6 * log10(inverse_wave_age);
else
  error('Inverse wave age outside of valid range of 0.84 < \Omega < 5');
end
BigGamma = exp( -( sqrt(k/kp) - 1).^2/(2*sigma^2));
Jp = littlegamma.^BigGamma;  % eqn (3)
LPM = exp( -(5/4) * (kp./k).^2);  %  eqn(2)


% Corrected to agree with eqn (32) - note LPM an JP are moved here.
% Fp =  LPM .* Jp .* exp( -(inverse_wave_age/sqrt(10)) *  ( sqrt(k/kp) - 1));

Fp = exp( -(inverse_wave_age/sqrt(10)) *  ( sqrt(k/kp) - 1));

Bl = 0.5 * alpha_p * (cp ./ c) .*  Fp;

%
% Short-wave curvature spectrum (sec. 5.2.2)
%

Fm = exp( -0.25 * ((k/km) - 1 ).^2);

if ( u_star < cm)
  alpha_m = (1e-2) * ( 1 + log(u_star/cm));
else
  alpha_m = (1e-2) * (1 + 3 * log(u_star/cm));
end

Bh = 0.5 * (alpha_m * cm ./ c ) .* Fm ;

%
% combine short and long wave for the omnidirectional spectrum.
%
S =  k.^(-3) .* Jp  .* LPM .* (Bh + Bl); % old version


%%%%%% WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%S =  k.^(-4) .*  Jp .* LPM .* (Bh + Bl); % fixed ?  version
%%%%% ACCORDING THE THEORY IT MUS BE -3.




return


