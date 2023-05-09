[x, y, z, clkdum] = ephem(eph, gpstime);
neph = size(eph,2);
xset(:,eph(1,:)) = x;
yset(:,eph(1,:)) = y;
zset(:,eph(1,:)) = z;
xgps = xset(:,satsvis(1,:));
ygps = yset(:,satsvis(1,:));
zgps = zset(:,satsvis(1,:));
for i=1:ntimes
  obsim = -pr_obs( [ 262013.0460 -4855185.5731  4114264.2173  0], ...
		   [xgps(i,:)' ygps(i,:)' zgps(i,:)'], 0)
  obsset(:,i) = obsim;
end
