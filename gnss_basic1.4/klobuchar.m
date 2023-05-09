function [iono_delay, zenith_delay, x] = ...
    klobuchar(azdeg, eldeg, gpstime, latdeg, londeg, alpha, beta)
%
% 11/04 - output the period x - so that I can verify some at night and
% some at day.
%
% Implements the Klobuchar "user" ionosphere correction defined for the
% GPS data message.
%
% Inputs el and az in DEGRESS, and the arrays of the alpha and beta
% coefficients. 
%
% Localtime is in seconds.  azdeg, eldeg and localtime must be the same
% dimensions
%
npts = size(gpstime,1);

% convert to semi-circles
AZ = azdeg / 180;
AZrad = azdeg * pi/180;
EL = eldeg / 180;
latu = latdeg / 180;
lonu = londeg / 180;

%
% Location of IPP
%
psi = 0.0137 ./ (EL + 0.11) - 0.022;

lati = latu + psi .* cos(AZrad);
lati(logical(lati > 0.416)) = 0.416;
lati(logical(lati < -0.416)) = -0.416;

longi = lonu + psi .* sin(AZrad) ./ cos( lati * pi);

%
% Convert to geomagnetic latitude
%
latm = lati + 0.064 * cos( ( longi - 1.617)*pi);

%
% local time (sec)
%
localtime = mod(4.32E4 * longi + gpstime, 86400);

%
% series expansions for amplitude and period
%
alpha = fliplr(alpha);
beta = fliplr(beta);
AMP = polyval(alpha, latm);
AMP = max(AMP, 0);
PER = polyval(beta, latm);
PER = max(PER, 0);
%
% Apply to half-cosine model
%

x = 2 * pi * (localtime - 50400)./PER;

dayflag = logical(abs(x) < 1.57);
zenith_delay = 5.0e-9 * ones(npts,1);
zenith_delay(dayflag) = zenith_delay(dayflag) + ...
    AMP(dayflag) .* ( 1 - x(dayflag).^2 /2 + x(dayflag).^4/24);

%
% Obliquity factor
%

F = 1.0 + 16.0 * (0.53 - EL).^3;

iono_delay = F .* zenith_delay;

return


