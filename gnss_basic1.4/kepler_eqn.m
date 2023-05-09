function [Fkep, Jkep] = kepler_equ(E, M, eccen);
%
% Kepler's equatio
%
Fkep = E - eccen .* sin(E) - M;
if (nargout > 1)
  Jkep = diag( 1-eccen.*cos(E));
end
return