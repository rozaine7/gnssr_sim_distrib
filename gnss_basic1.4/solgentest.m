%
% Test script for solngen.m
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
Solnpt = 1;
PV0 = [x0, y0, z0, t0];
%PV0 = PVstate;
%%gpsxyztest = squeeze(gpsxyz(Solnpt,:,:))';
%%prtest = prc(Solnpt,:)';

%%rawprtest = pr(Solnpt,:)';

%%goodpr = logical( rawprtest ~= -999);
%%prtest = prtest(goodpr);
%%gpsxyztest = gpsxyztest(goodpr,:);

%%preres = pr_obs(PV0, gpsxyztest, prtest);

PVstate = PV0;
for k=1:10
   res = pr_obs(PVstate, gpsxyztest, prtest)
   G = gmatrix(PVstate(1:3), gpsxyztest);
   dstate = inv(G' * G) * (G') * res
   PVstate = PVstate + dstate'
end
%options = optimset('TolX',1e-15,'TolFun',1e-15,'Jacobian','on');
%PVstate = lsqnonlin('pr_obs',PV0,[],[],[], gpsxyztest,prtest);

postres = pr_obs(PVstate, gpsxyztest, prtest);