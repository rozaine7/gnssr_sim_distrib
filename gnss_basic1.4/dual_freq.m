function [code_delay, carrier_delay] = dual_freq( rhoL1, rhoL2, phiL1, phiL2)
%
%  Caulculation of code and carrier ionospheric delays (at L1) 
%  from measurements made at L1 and L2.
%
f1 = 1575.42e6;
f2 = 1227.60e6;
c = 2.99792458E8;
lambda1 = c/f1;
lambda2 = c/f2;

code_delay = (f2^2)/(f2^2 - f1^2) * (rhoL1 - rhoL2);
carrier_delay = (f2^2)/(f2^2 - f1^2) * (lambda2*phiL2 - lambda1*phiL1);

return