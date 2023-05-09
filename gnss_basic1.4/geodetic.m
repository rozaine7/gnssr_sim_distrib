function [latg, longg,  htg] = geodetic(X, Y, Z);
%
%  iterative solution for geodetic latitude
%
a = 6378137;
b = 6356752.314;
eccen2 = (a^2 - b^2)/a^2
p = sqrt(X^2 + Y^2);
htg0 = 0;
latg1 = pi;  % something large to force at least one iteration
latg0 = atan( (Z/p) / ( 1- eccen2));
while abs(latg1 - latg0)  > 1e-10
    N = a^2 / sqrt( a^2 * cos(latg0)^2 + b^2 * sin(latg0)^2 );
    htg = p / cos(latg0) - N;
    latg = atan( (Z/p) * ( 1 - eccen2 * N / ( N + htg))^(-1) );
    latg1 = latg0;
    latg0 = latg;
end

longg = atan2 ( Y/p, X/p);

return
  
