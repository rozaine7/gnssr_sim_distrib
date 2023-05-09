function [GPSWK, GPSD] = YMD2GPSW(YYYY, M, D, span);
%
% [GPSWK, GPSD] = YMD2GPSW(YYYY, M, D, span)
% 
% Computes GPS Week (GPSWK) and day in the week (GPSD) 
% given 4-digit year (YYY), Month (M) and Day (D).  
% To accomodate batch processing, the outputs will be matrices computes for 
% each day from D-span to D+span. (span optional)
%
  if (nargin < 4)
    span = 0;
  end
MJD0 = YMD2MJD(1980, 1, 6, 0);
day_range = [-span:0 1:span];
MJD_mid = YMD2MJD(YYYY, M, D, 0);
MJD = MJD_mid + day_range';

GPSWK = floor( (MJD - MJD0)/7);
GPSD = mod( MJD - MJD0, 7);

return
