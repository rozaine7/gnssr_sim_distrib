function [tau, resnorm, residual] = chartime( lags, autocorr)
%
% Compute the correlation time (characteristic time) from a Gaussian model
% fit to an autocorrelation function
%
ndelays = size(autocorr,2);
tau = zeros(ndelays,1);
var = zeros(ndelays,1);
resnorm = zeros(ndelays,1);
residual = zeros(ndelays,1);

for kdelay = 1:ndelays
  [xest, resnorm(kdelay), res] = lsqnonlin('gauss_corr',[1, 1], [], [],[],autocorr(:,kdelay), lags');
  residual(kdelay) = sum(res);
  var(kdelay) = xest(1);
  tau(kdelay) = xest(2);
end

end




