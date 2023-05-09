function [XE, YE, ZE] = ecef_orb_rot( X, Y, Z, timelag)
%
% Rotates the GPS satellite positions from an ECEF reference at time 
% t - timelag to an ECEF reference at time t, so that the range can be
% computed in the same reference as the receiver.
%
c = 2.99792458e8;
omegaE = 7.2921151467e-5;  

cw = cos(omegaE * timelag);
sw = sin(omegaE * timelag);

XE = cw * X + sw * Y;
YE = -sw * X + cw * Y;
ZE = Z;

return
