function [lat,long,alt] = WGS84_2_latlong(X, Y, Z)
%
%  this function converts WGS84 ECEF coodinates to lat-longs
%
%  using the algorithm from Hoffman-Wellenhopf p 232
%
a = 6378137;
b = 6356752.314;
eprime2 = (a^2-b^2)/b^2;
e2 = (a^2-b^2)/a^2;
p = sqrt(X.^2+Y.^2);
theta = atan( (Z*a)./(p*b));
lat = atan((Z+eprime2*b*sin(theta).^3)./...
      (p-e2*a*cos(theta).^3) );
long = atan2( Y./p, X./p);
N = a^2 ./ sqrt(a^2 .* cos(lat).^2 + b^2 .* sin(lat).^2 );
alt = p./cos(lat) - N;
return
