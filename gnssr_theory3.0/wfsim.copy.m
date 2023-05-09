function [pcdarray, y, ftilde, corrspec, fdsurf, Rtau, Rtau_fix, yindep] = wfsim(mp, cp, CN0, nwf, nic, binbinflag);
%
% GNSS-R waveform generator 
%

% 
% Compute the spectra for the time series from each bin.
%
fprintf(' Computing the power spectra ... \n');

dftilde = mp.dftilde; 
ftilde = [-mp.ftilde_max:dftilde:mp.ftilde_max];  % Hz
nftilde = size(ftilde,2);


for k=1:size(ftilde,2)
 [wf_array, pcdarray, corrspec(:,k), pwr, fs, fdsurf] = ...
     wf_from_gcp_dopp( 0, cp, mp, 0, ftilde(k), dftilde);
 fprintf('   ftilde = %12.5f Hz \n', ftilde(k))
end

if (binbinflag ~=0) 
       [wf_array, pcdarray, dummy, pwr, fs, fdsurf, Rtau] = ...
          wf_from_gcp_dopp( 0, cp, mp, 0, 0, dftilde);

      size(wf_array)
end

yindep = datasim(ftilde', corrspec', cp.Ti, nwf, nic, CN0);

if (binbinflag ~=0)
%
% Bin-bin correlation
%

%
% This is a bit of a hack
%
   Rtau_fix = covfix(Rtau, 300, 5e-3, 1);

   ynorm = yindep./(ones(1000,1) * wf_array');

   y=ynorm * Rtau_fix;
else

   y = yindep;
   
end

return
