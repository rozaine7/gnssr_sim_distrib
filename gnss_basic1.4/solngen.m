function [rxyz, tr] = solngen ( gpsxyz_all, pr_all, rxyz0)
%
% Generate GPS receiver point position using pseudrange measurements
% and GPS satellite positions.
%
% pr and gpsxyz are both equi-dimensioned arrays. Satellites can be
% selected and the pr can be smoothed prior to calling this subroutine.
%
% gpsxyz is the GPS positions computed at the time of signal
% transmission, in the ECEF frame at the instant the signal is received.  
% ie., this rotation is performed prior to passing the ECEF satellite
% positons to this subroutine.
%
% It is assumed that all corrections (clock biases, iono, etc) are
% applied to the pr array prior to calling this subroutine.
%
[ntimes, three, nsv] = size(gpsxyz_all);
rxyz = zeros(ntimes,3);
tr = zeros(ntimes,1);
for i=50:50
  PV0 = [rxyz0, 0];
  gpsxyz = squeeze(gpsxyz_all(i,:,:))';
  pr = pr_all(i,:)';
  size(gpsxyz)
  size(pr)
  goodpr = logical( pr > 0);
  pr = pr(goodpr);
  gpsxyz = gpsxyz(goodpr,:);
  pr
  gpsxyz
  PVstate = lsqnonlin('pr_obs',PV0,[],[], ...
            optimset('TolX',1e-15,'TolFun',1e-12),gpsxyz,pr);
  rxyz(i,:) = PVstate(1:3);
  tr(i,1) = PVstate(4);
end