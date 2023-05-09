function [o_minus_c, J] = pr_obs(rcvr_PT, sat_xyz, obsset)
%
% GPS pseudorange measurement equations.  Takes a vector set of 
% Nobs PR observations (obsset) in, along with Nobs X 4 set of satellite
% ECEF positions.  Outputs observed-computed pr and the Jacobian (as
% req'd for the MATLAB fsolve function).
%
nsat = size(sat_xyz,1);
o_minus_c = obsset - (sqrt((rcvr_PT(1)-sat_xyz(:,1)).^2 + ...
		          (rcvr_PT(2)-sat_xyz(:,2)).^2 + ...
		          (rcvr_PT(3)-sat_xyz(:,3)).^2)+rcvr_PT(4));
if (nargout > 1)
  los = sat_xyz - ones(nsat,1) * rcvr_PT(1:3);
  range = sqrt( sum(los.^2,2))*ones(1,3);

  J = [los ./ range, -1*ones(nsat,1)];
end

return
						  
