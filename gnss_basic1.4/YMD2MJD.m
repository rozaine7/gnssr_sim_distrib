function MJD = YMD2MJD (Y, M, D, H)

% this function converts buoy calender time into Modified Julian Date

before_march_1 = logical (M <= 2);
ysmall = Y;
msmall = M;
ysmall(before_march_1) = Y(before_march_1) - 1;
msmall(before_march_1) = M(before_march_1) + 12;
B = floor(ysmall/400) - floor(ysmall/100) + floor(ysmall/4);
MJD = 365*ysmall-679004+B+floor(30.6001*(msmall+1))+D + H /24;   

return