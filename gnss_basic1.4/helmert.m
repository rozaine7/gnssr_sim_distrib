function x2 = helmert( x1, c, mu, alpha1, alpha2, alpha3)
%
% Helmert transformation for cartesian rotations
%
T1 = rotmat(1, alpha1);
T2 = rotmat(2, alpha2);
T3 = rotmat(3, alpha3);
x2 = c + mu * T3 * T2 * T1 * x1;

return