function [E, N, U, ECEF_2_ENU_mat] = ECEF_2_ENU ( X, Y, Z, phi, lambda)
%
% Rotates ECEF cartesian position to ENU position.
%
% X, Y, Z can be matricies, phi, lambda in radians.
%
ECEF_2_ENU_mat = ...
  [ -sin(lambda), cos(lambda), 0; ...
  -sin(phi)*cos(lambda), -sin(phi)*sin(lambda), cos(phi); ...
    cos(phi)*cos(lambda), cos(phi)*sin(lambda), sin(phi)];
E = -sin(lambda)*X + cos(lambda)*Y;
N = -sin(phi)*cos(lambda) * X - sin(phi)*sin(lambda) * Y + cos(phi)*Z;
U = cos(phi)*cos(lambda) * X + cos(phi)*sin(lambda) * Y + sin(phi)*Z;

return
