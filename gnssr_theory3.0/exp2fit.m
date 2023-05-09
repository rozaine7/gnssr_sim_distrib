function [Serr, Jac] = exp2fit ( params, f, x )
%
% Residual function for fitting Ae^-(x/B)^2 function shape to power
% spectrum
%
A = params(1);
B = params(2);

Serr = x - A * exp(-(f/B).^2);

Jac = [ -exp(-(f/B).^2), -(A/B^2)*exp(-(f/B).^2)];

end

