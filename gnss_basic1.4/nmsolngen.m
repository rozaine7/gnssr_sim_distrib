[nsats, ntimes] = size(obsset);
[x, y, z, clkdum] = ephem(eph, gpstime);
neph = size(eph,2);
xset(:,eph(1,:)) = x;
yset(:,eph(1,:)) = y;
zset(:,eph(1,:)) = z;
xgps = xset(:,satsvis(1,:));
ygps = yset(:,satsvis(1,:));
zgps = zset(:,satsvis(1,:));
PV = zeros(ntimes,4);
dx = sv_clock_fix(eph, sincebom);
dxset(:,eph(1,:)) = dx;
dxgps = dxset(:, satsvis(1,:));
for i=1:1
%  obs = obsset(:,i)-dxgps(i,:)';
  obs = obsset(:,i);
  gpsxyz = [xgps(i,:)', ygps(i,:)', zgps(i,:)'];
  PVstate = lsqnonlin('pr_obs',PV0,[],[], ...
            optimset('TolX',1e-15,'TolFun',1e-12),gpsxyz,obs);
  PV(i,:) = PVstate;
end