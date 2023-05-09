function [x, y, r, theta, x0, xs, a, b] = isorange( alt, gammadeg, delay)
%
% compute bistatic isorange lines assuming a flat Earth.
%
% This uses the full expression from Beckmann & Spizzochino without any
% approximation.  x is cartesian coordinate along incidence plane, y is
% perp. to incidence plane. The origin is the "apparent" receiver nadir
% point over a flat Earth.  r, theta are polar coordinate used to
% generate this locus.  x0 is location of center of isorange ellipse
% xs is specular point.
%
% Inputs:  gamma and alt are incidence angle (deg, measured from specular
% point), alt is receiver altitude (m). delay is path delay relative to
% the specular point (in meters)
%
% when used with a set of delays, they are expected to be in a row
% vector. Columns of outputs all correspond to the same delay.
%
% Angle theta measured from the specular point - this is the *Same* as in
% the integrand: dp_dtheta_dopp_sumf.c.  As of 2/15 - it is NOT the same as
% in isodopp
% 
gamma = gammadeg * pi/180;

xs = alt/tan(gamma);
x0 = (cos(gamma)/sin(gamma)^2) * delay + xs;

theta = [-pi:0.01:pi]' * ones(1,size(delay,2));

delta_prime = delay + alt * sin(gamma);

C = sqrt( ones(size(theta,1),1) * (delta_prime/sin(gamma)).^2 - alt^2);


r = C ./ sqrt( sin(theta).^2 + (sin(gamma)^2) * cos(theta).^2 );

x = ones(size(theta,1),1) * x0 + r .* cos(theta);
y = r .* sin(theta);

%
% Semimajor and semi-moinor axes of ellipse
%

a = sqrt(2 * delay * alt * sin(gamma))/sin(gamma)^2;
b = sqrt(2 * delay * alt * sin(gamma))/sin(gamma);

return