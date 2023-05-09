function rt = td_corr(code1, code2, istart, lag, nsamples)
%
% Time-domain calculation of the autocorrelation between two code
% sequences.  This applies the following definition exactly:
%
% rt(lag) = (1/N) * \sum_0^N code1(i) * code2(i+lag)
%
% This is not efficient, but is good for demonstrating this basic
% principle.
%
% This also assumes that the lengths of code1 and code2 are at least
% istart + nsamples + lag.
%
% code1 and code2 must be column vectors, lag must be a row vector.
% rt will be a row vector, with each element corresponding to a value of lag
%
istart
nsamples
nlags = size(lag,2)
code2_lag = zeros(nsamples, nlags);
for i=1:nlags  % there is probably a more clever way of doing this.
  code2_lag(:,i) = code2(istart+lag(i):istart+nsamples-1+lag(i));
end
rt = sum( code1(istart:istart+nsamples-1)*ones(1,nlags) .* ...
	  code2_lag,1)/nsamples;

return
      
