function pr = pr_obs(rcvr_PT, sat_xyz, obsset)
%
% GPS pseudorange measurement equations.  Takes a vector set of 
% Nobs PR observations (obsset) in, along with Nobs X 4 set of satellite
% ECEF positions.  Outputs observed-computed pr and the Jacobian (as
% req'd for the MATLAB fsolve function).
%
nsat = size(sat_xyz,1);
pr = obsset - (sqrt((rcvr_PT(1)-sat_xyz(:,1)).^2 + ...
                    (rcvr_PT(2)-sat_xyz(:,2)).^2 + ...
                    (rcvr_PT(3)-sat_xyz(:,3)).^2)+rcvr_PT(4));

return
						  
