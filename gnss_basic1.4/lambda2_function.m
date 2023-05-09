function [lambda, lambda2] = lambda2_function(tau_chips);
%
%  computes the lambda^2 function in units of half-code chips
%
npts = size(tau_chips,1);
lambda = max( ones(npts,1)-abs(tau_chips), zeros(npts,1));
lambda2 = lambda.^2;

return