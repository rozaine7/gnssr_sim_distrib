function [xecef, yecef, zecef, sv_clock, gpstime] = alman(alm, gpstime)
%
%  Program to read in the GPS satellite almanac data and
%  compute a cartesian position from it.
%
%  Input:  alman = 12 X Nsat = broadcast almanac for N satellites.
%          gpstime = Ntime X 1 - time desired for almanac is GPS time
%          sec.
%  Output:  xecef, yecef, zecef, sv_clock - ea. Ntime X Nsat - ECEF
%  positons and clock offsets of GPS sats at times in gpstime.
%
%  Calls the ephem.m code, with some of the parameters set to zero
%

nsat = size(alm,2);
eph = zeros(23,nsat);

eph(1,:) = alm(1,:);
eph(2,:) = alm(3,:);
eph(3,:) = alm(4,:);
eph(4,:) = alm(5,:);
eph(5,:) = alm(6,:);
eph(6,:) = alm(7,:);
eph(7,:) = alm(8,:);
eph(8,:) = alm(9,:);
eph(9,:) = alm(10,:);
eph(21,:) = alm(11,:);
eph(22,:) = alm(12,:);
%
% Defined Ephemeris data are read from input matrix
%


[xecef, yecef, zecef, sv_clock] = ephem(eph, gpstime);

return