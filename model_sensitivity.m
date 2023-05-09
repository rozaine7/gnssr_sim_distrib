function [ dtbvdmss, dtbhdmss ] = poly4model( mss )
% 
% 4th order polynomial fit - found by JLG with better agreement than the
%  piecewise fit JV did.
%
% Coefficients found from MATLAB built in POLYFIT function
%
dtbvdmss =  -2.6185e+08 * 4 * mss.^3 + 2.7042e+07 * 3 * mss.^2 -1.0072e+06 * 2 * mss + ...
 1.6221e+04; 

dtbhdmss = -2.8409e+08 * 4 * mss.^3 + 2.9142e+07 * 3 * mss.^2 -1.0782e+06 * 2 * mss + ...
 1.7310e+04; 

%
% Old version - in terms of sigmau^2 = sigmac^2 = sigma_iso^2
%
%tbv = -4.1896e+09 * mss.^4 + 2.1634e+08 * mss.^3 -4.0288e+06 * mss.^2 + ...
% 3.2442e+04 * mss + 25.2696; 

%tbh = -4.5454e+09 *mss.^4 + 2.3314e+08 * mss.^3-4.3127e+06 * mss.^2 + ...
% 3.4620e+04 * mss -27.5029; 
end

