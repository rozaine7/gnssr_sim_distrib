function [wfpwr, chip_delays, Doppler_freqs, pwr_array, fs, dp_real, dp_imag] = ...
    d2map ( case_params, model_params );
%
%  1/04 modified to use structures to pass all of the parameters through
%  case_params - anything specific to one scenario
%  model_params - defining the type of model used, and necessary
%  constants
%
%  This function calls the waveform integration over a range of 
%  Doppler bins to generate a 2D waveform
%
%  Requires ACF to be generated first, and saved in the global 
%  variable acf2matrix 
%
%  Doppler is being re-introduced in the MEX-file version - note that
%  NaN's result when running with VR=VG=0 identically.  model_flag won't 
%  be used to turn Doppler off - to simplify the code.  For zero-Doppler
%  examples, will have to set the velocity very small.
%
% wfpwr, chip_delays, Doppler_freqs, pwr_array, fs] = ...
%    d2map ( case_params, model_params );
%

Ti = case_params.Ti;
alt_m = case_params.alt_m;
chipsize   = model_params.acfmodel.chipsize;
nDoppbins = case_params.nDoppbins;
maxdelay  = case_params.maxdelay;
gamma = case_params.gammadeg * pi/180;
doppbin_frac = model_params.doppbin_frac;

deltaf = doppbin_frac / (2*Ti); % spacing between dopp bins 
Dopp_freq_array = [-nDoppbins*deltaf:deltaf:deltaf*nDoppbins];
nDopp = size(Dopp_freq_array,2);

%
% Removed ACF generation from here - replaced it with lookup of
% global array
%

%   alt = alt_m; not sure why this is needed.


for i=1:nDopp
  [wf_array, pcdarray, cspec, pwr, fs, dmy1, dmy2, dmy3, dmy4, dp_real, dp_imag] = wf_from_gcp_dopp ...
       ( Dopp_freq_array(i), case_params, model_params, 0, 0, 0);
  wfmask = logical( (pcdarray >= - 1.5 * chipsize) & ...
                    (pcdarray < maxdelay) );
  chip_delays(:, i) = pcdarray(wfmask)';
  wfpwr(:,i) = wf_array(wfmask);
  pwr_array(:,i) = pwr;
end

nwfsamples = size(wf_array,1);

Doppler_freqs = ones(nwfsamples, 1) * Dopp_freq_array;

return







