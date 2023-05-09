function [xyzstate, G, dx] = point_soln(xyzstate0, gpsxyz, obsset)
%
%  Iterative least squares GPS point solumn
%
oc = pr_obs(xyzstate0, gpsxyz, obsset);
G = gmatrix(xyzstate0(1:3),gpsxyz);
dx = (inv(G'*G) * G') * oc;
xyzstate = xyzstate0' + dx;

return