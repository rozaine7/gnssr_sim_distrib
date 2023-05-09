function [pcdarray, y, ftilde, corrspec, fdsurf, Rtau_fix, wf_array, Rtau0, ...
    comptime, simtime, compcpu, simcpu, ntheta] = ...
    wfsim(mp, cp, CN0, nwf, nic, binbinflag, Bmin, powervariation)
%
% Synthtic GNSS-R waveform generator 
%  
% Output:
%   pcdarray = delay (m) [1 x no. of bins]
%   y = waveform (non-dimensiona1).   
%      If nic = 0 - this will nwf realizations of the complex (I,Q) from coherent integratio
%        of cp.Ti sec. 
%      If nic > 0 - this will nic realizations of the power (I^2+Q^2) from the incoherent
%        sum of nic waveforms, each with coherent integration of cp.Ti
%   ftilde = (Hz) [1 X nftilde] frequency axis for spectral model of delay-bin time
%        series.  nftilde set by mp.ftilde_max and mp.dftilde
%   corrpsec = (1/Hz) [ no. bins X nftilde] Power spectrum of each
%        delay-bin time series, computed at the frequencies given in
%        ftilde.
%   fdsurf = presently unused
%   Rtau_fix = [no.bins X no.bins] non-dimensional - covariance matrix of
%        the power vs. delay waveform. 
%   wf_array = (non-dimensiponal) - numerical model for power vs. delay
%        waveform, at delays given in pcdarray.
%
% Input:
%   mp = "model parameters" structure passed to the scattering model to
%     define the basic model assumptions.  See "wf_from_gcp_dopp.m" for
%     list of elements
%   cp = "case parameters" structure passed to the scattering model to
%     setup a specific problem or case.   See "wf_from_gcp_dopp.m" for
%     list of elements
%   CN0 = Carrier to noise ratio (Hz) of incident signal (at present, the
%     surface reflectivity is not included in the scattering model, so this
%     must be included first by reducing CN0.  This should be fixed.
%   nwf = number of waveform realizations to be generated, sequentially in
%     time, each with coherent integration time of cp.Ti. If nic=0, the
%     complete set of these complex waveforms are returned in the matrix y.
%     If nic >0, then they are incoherently summed, to produce one of the
%     nic power waveforms, this is repeated nic times. 
%   nic = number of power waveforms generated, each formed by the
%     incoherent sum of nwf voltage waveforms. If nic =0, the waveforms are
%     not incoherently averaged, but returned as an array of complex
%     voltage waveforms. 
%   binbinflag = 1 - generate waveforms that have the correct correlation
%     between delay bins. 
%              = 0 - generate independent time series for each delay bin,
%              = -1 - compute bin-bin covariance matrix, but do not
%              normalize it - this is just to compute the matrix, not
%              simualtion, as the signal power in each bin is set by the
%              integrated spectrum.  
%   Bmin- minimum bandwidth of signal due to the de-correlation of surface
%     (Hz).  Setting Bmin=0 will not modify the spectrum.  
%   powervariation = Arbitrary gain can be implemented to simualate 
%     calibration errors in total power          
%
% Compute the spectra for the time series from each bin.
%
fprintf(' Computing the power spectra ... \n');

%% Begin timing of model computation

%
% set equal to zero, in case you dont run the simulation
%

simcpu = 0;
simtime = 0;
y =-999;
Rtau_fix = -999;

tic  
compstart = cputime;

dftilde = mp.dftilde; 
ftilde = [-mp.ftilde_max:dftilde:mp.ftilde_max];  % Hz
nftilde = size(ftilde,2);


for k=1:size(ftilde,2)
 [wf_array, pcdarray, corrspec(:,k), pwr, fs, fdsurf] = ...
     wf_from_gcp_dopp( 0, cp, mp, 0, ftilde(k), dftilde);
 fprintf('   ftilde = %12.5f Hz \n', ftilde(k))
end

[corrspec, Bfix, B, S0] = specfix(ftilde, corrspec, Bmin);


if (binbinflag ~=0) 
  fprintf('Computing Bin-bin correlation matrix ... \n')
  %
  % note - do not put the ntheta argument for the previosu one - otherwise
  % you will generate a covaraince matrix at each point in the spectrum!
  %
  [wf_array, pcdarray, dummy, pwr, fs, fdsurf, Rtau,Rtau0, ntheta] = ...
          wf_from_gcp_dopp( 0, cp, mp, 0, 0, dftilde);
      

%
% Bin-bin correlation
%
%  Normalize spectrum to have a area (ie. total power) of unity.  This will
%  preserve the correlation time of the signal, when the correlation
%  between the individual delay bins is introduced by multiplying by the
%  Rtau matrix. 
%
%  Necessary to remove cases in which the integration of the spectrum is
%  zero (which may be do to finite length of ACF - still investigating
%  this) - otherwise you will get NaN's when the synthetic data is
%  generated
%
   
   spec_int = trapz(ftilde, corrspec,2);
   spec_int(logical(spec_int == 0)) = 1; % Added to eliminate NaN's
                                         % under these cases, the denom. does 
                                         % not matter as the power is zero.
   
   corrspec_unity = corrspec ./ (spec_int * ones(1,size(ftilde,2)));
%% End timing of model computation 
  comptime =  toc;
  compcpu = cputime-compstart;
   
%% Begin timing of synthetic data generation
tic

if nwf ~= 0 
simstart  = cputime;

   if (binbinflag ~= -1)  % -1 is test value to not normalize Rtau 
      Rtau_fix = covfix(Rtau, 300, 1e-5, 1e-10, 1);
      fprintf('Generating synthetic waveforms ... \n')
     % save test
      y = datasim(ftilde', corrspec_unity', cp.Ti, nwf, nic, CN0, powervariation, Rtau_fix, Rtau0);
   else
      fprintf(' Generating non-normalize R(tau,tau) - DO NOT USE FOR SYNTEHTIC DATA \n')
      Rtau_fix = Rtau;
      y = 0;
   end
%%    End timing of synthetic data gneration
simtime = toc
simcpu = cputime - simstart;
end

else
    
    
%
% Each delay bin is independent
if nwf ~= 0 
   fprintf('Generating synthetic waveforms, bin-bin uncorrelated  ... \n')
   y = datasim(ftilde', corrspec', cp.Ti, nwf, nic, CN0, powervariation);
   Rtau_fix = eye(size(y,2));
%   Rtau = 0;

end
end


return
