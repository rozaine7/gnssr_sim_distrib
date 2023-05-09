function [sigma2_u, sigma2_c, sigma2omni] = sigma_from_spec( ...
                              U10mat, inverse_wave_age, kmaxmat)
%
% 4/04 - allow column vector for U10 and outputs.
% Integrates wave spectrum up to wave number kmax to give the mss.
%
npts = size(U10mat,1);
if size(kmaxmat,1) == 1
  kmaxmat = kmaxmat * ones(npts,1);
end
sigma2_u = zeros(npts,1);
sigma2_c = zeros(npts,1);
for ipt = 1:npts
  U10 = U10mat(ipt);
  kmax = kmaxmat(ipt);
  dSx = inline( ...
    '(cos(phi)^2) * k.^3 .* elfouhaily2d( k, phi, U10, inverse_wave_age)', ...
    'k','phi','U10','inverse_wave_age');
  dSy = inline( ...
    '(sin(phi)^2) * k.^3 .* elfouhaily2d( k, phi, U10, inverse_wave_age)', ...
    'k','phi','U10','inverse_wave_age');
  sigma2_u(ipt,1) = dblquad(dSx, 1e-3, kmax, -pi, pi, ...
		  1e-8,[], U10, inverse_wave_age);
  sigma2_c(ipt,1) = dblquad(dSy, 1e-3, kmax, -pi, pi, ...
		  1e-8,[], U10, inverse_wave_age);
  dSo = inline( ' k.^2 .* elfouhaily1d( k, U10, inverse_wave_age)', ...
    'k','U10','inverse_wave_age');
  sigma2omni(ipt,1) = quad(dSo, 1e-3, kmax, [],[], U10, inverse_wave_age);
end
return