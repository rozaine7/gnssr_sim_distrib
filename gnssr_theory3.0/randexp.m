function p = randexp ( pbar, ncorr, nelements );
%
%  Generates an exponentially distributed random variable 
%  with mean and std. of pbar and a correlation time of ncorr.
%
%  This is done by:
%
%   1.) Generating two gaussian RV's  with zero mean and unity std. 
%   2.) Passing them through a moving average filter with window ncorr
%   3.) Scaling theses results by amplitde * sqrt(ncorr);
%   4.) Computing the exponentially distributed power as the sum of squares of %       these two RV's.
%
nsamples = size(pbar,1);
I0 = randn(nsamples, nelements+ncorr);
Q0 = randn(nsamples, nelements+ncorr);
If = filter([1/ncorr], [1, 1/ncorr-1], I0, [ ], 2);
Qf = filter([1/ncorr], [1, 1/ncorr-1], Q0, [ ], 2);
I = If(:,floor(ncorr/2)+1:floor(ncorr/2)+nelements);
Q = Qf(:,floor(ncorr/2)+1:floor(ncorr/2)+nelements);
%mean_pwr = mean((I.^2 + Q.^2), 2) * ones(1, nelements)
p = (pbar * ones(1,nelements)) .* ((I.^2 + Q.^2)/2);
return

