function [X, Y, Z, N] = latlong_2_WGS84(lat, long, alt)
%
%  converts from latitude, longitude and altitude to cartesian
%  position in the WGS-84 system
%
a = 6378137;
b = 6356752.314;
N = a^2 ./ sqrt(a^2 .* cos(lat).^2 + b^2 .* sin(lat).^2 );
X = (N + alt) .* cos(lat) .* cos(long);
Y = (N + alt) .* cos(lat) .* sin(long);
Z = ( b^2 * N /a^2 + alt).*sin(lat);
return

