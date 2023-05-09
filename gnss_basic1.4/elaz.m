function [azd, eld] = elazcomp( xt, yt, zt, x, y, z)
  size(xt)
  size(x)
%
%  compute the elevation and azimuth of a satellites.
%
%  xt, yt, zt are nX1
%
%  x, y, z are nX31 array of satellite WGS84 positions.
%
%  azd, eld  will be n X 31
%
npts = size(x,1);
nsv = size(x,2);
az = NaN * ones(npts,31);
el = NaN * ones(npts,31);
   
radt = sqrt(xt.^2 + yt.^2 + zt.^2);
north = [ -xt.*zt, -yt.*zt, xt.^2+yt.^2];
radn = sqrt( dot(north, north));

xtm = xt*ones(1,nsv);
ytm = yt*ones(1,nsv);
ztm = zt*ones(1,nsv);


[latr, lonr, htr] = WGS84_2_latlong( xt, yt, zt);
latr = latr * ones(1,nsv);
lonr = lonr * ones(1,nsv);
htr = htr * ones(1,nsv);

dxGPS =  x-xtm;
dyGPS = y-ytm;
dzGPS = z-ztm;

OTH = logical(dxGPS.*xtm + dyGPS .* ytm + dzGPS .*ztm > 0);
radgps = sqrt( dxGPS.^2+dyGPS.^2 + dzGPS.^2);
             
%el(OTH) = asin( ...
%(xt*dxGPS(OTH)+yt*dyGPS(OTH)+zt*dzGPS(OTH))./(radgps(OTH)*radt) );



%
% Better way - by using rotations - first to the primes - which are in
% the equatorial plane
%

dxprime = dxGPS.*cos(lonr)+dyGPS.*sin(lonr);
dyprime = -dxGPS.*sin(lonr) + dyGPS.*cos(lonr);
dzprime = dzGPS;

% 
% Now to the N,U,W system
%

N = -dxprime.*sin(latr) + dzprime.*cos(latr);
U =  dxprime.*cos(latr)+dzprime.*sin(latr);
E = dyprime;

sinel = U ./ radgps;

el(OTH) = asin(sinel(OTH));

proj_horiz = sqrt( N.^2 + E.^2);

sinaz = E ./ proj_horiz;
cosaz = N ./ proj_horiz;


az(OTH) = atan2( sinaz(OTH), cosaz(OTH));

eld = el * 180/pi;
azd = az * 180/pi;

negaz = logical( azd < 0);
azd(negaz) = azd(negaz) + 360;

return
