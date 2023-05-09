%
% Script to generate map of ambiguity function
%
f = -5000:100:5000;
nfreq = size(f,2);
corrmat = zeros(5715,nfreq);

for k = 1:nfreq
  [lags, corrpwr] = ...
      code_search(5.714286e6, 1.405e6+f(k), hw2_data, pncodes, prn);
  corrmat(:,k) = corrpwr;
  k
end
				