function [yy, post_corr_delay, corr_spec, pwr, fs, fdopp, Rtau, Rtau0, ntheta, dp_real, dp_imag]  = ...
     waveform( doppoff, case_params, model_params, ttilde, ftilde, dftilde);
%
% 1/15 - modified to write up surface doppler vs. theta for plots in the
% paper. 
% 1/13 - correction of conv - replacement with explicit fft/ifft
% 10/12 - added comments to distribute to CYGNSS team.
% 6/10 - back to this to fix it for the paper and the PO proposal
%   not sure if bin-bin correlation is working yet.
% 6/08 Latest modifications, to generate bin-to-bin covariance matrix.
% 1/04 modified to use structures to pass parameters
%
% Master code to generate model waveforms, their correlation and spectra in
%    time and the bin-bin covariance matrix.
% Generate the reflected GPS waveform, given the constants of the
% Gram-Charlier PDF series.
% 
% Generates time- and frequency domain numerical models for:
%  R(tau, f, ttilde) = < Y(tau, f, tilde) Y^*(tau, f, t+ttilde) > 
%  S(tau, f, ftilde) = \int R(tau, f, ttilde) exp( -j*ttilde*ftilde) d ttilde
%
% In which Y is the complex correlation of the received signal x(t) over
%   coherent integration time Ti
%   Y(tau, f, ttilde) = \int_Ti x(t+tau) s(t) exp( -j*t*f) dt
%
% s(t) is the baseband local signal generated by the PRN code. 
%
% Output:
%   yy = R(tau, f, ttilde), when ttilde=0 this reduces to the waveform model 
%      I^2 + Q^2 which should be real.  
%   post_corr_delay = tau (m) delay axis for waveform
%   corr_spec = S(tau, f, ftilde) = correlation spectra of delay/doppler bin (tau,f).
%   pwr = intermediate result in computing yy - single integral over azimuthal 
%     coordinate - when convolved with auto-correlation to produce waveform 
%     (ttilde = 0) or autocorrelation
%     (ttilde ~= 0)
%   fs = Doppler shift at the specular point (Hz)
%   corr_spec_delta = not used in normal operation - test ouptut of the
%     numerically-implemented dirac delta function in the integration
%   Rtau = waveform covariance matrix - this will not be generated unless
%     it is given as an output - this is computationally intense. 
%   Rtau0 = covariane matrix for theremal noise, due to correlation of same
%     process with multiple delays. 
% 
% Inputs:
%   dopoff = doppler offset, will compute the Doppler at f+dopoff.  Not
%      presently used
%   case_params = "case parameters" structure - used to pass numbers for a specific 
%       case (scattering geometry, surface conditions, etc) to the model
%       wdirdeg = wind direction (deg) measured in right-hand sense from
%          scattering plane
%       alt_m = receiver altitude (m) from tangent plane
%       gammadeg = elevation of satellite (deg) measured from tangent plane
%       Ti = coherent integration time (sec) 
%       nDoppbins = number of Doppler bins to generate. Total doppler
%          dimension will be 2*nDoppbins + 1, using equally-space +/-
%          Dopplers, in increments of 1/(2*TI)
%       maxdelay = maximum delay to generate waveform (m)
%       prn = satellite PRN - allows use of full C/A code model produced by
%          acfgen. 
%       VR = receiver velocity (m/s) in scattering plane coordinates. 
%       VG = GPS satellite velocity (m/s) in scattering plane coordinates
%       PDF_params = Coefficients of the Gram-charlier model for slope PDF
%          [upwind variance, cross-wind variance, X X X X X]
%       roll = aircraft attitude- needed for accurate implementation of
%          antenna gain
%       pitch
%       yaw
% model_params = "model parameters" structure - all of the numerical values to
%   define the assumptions in the model. 
%       dftilde = increment in ftilde in generating spectrum
%       ftilde_max = range in ftilda (-ftilde_max to ftilde_max)
%       ttilde_max = range in generating the time-correlation of each
%          delay-bin time series (-ttilde_max
%       acfmodel = structure defining parameters of the code autocorrlation
%          model - see acfgen.m for details on what is in this
%       fullmodel = structure to pass parameters defining the scattering
%          model.  At present only one fullmodel.PDF, model for PDf of surface
%       wavelength = carrier wavelenght (m) 
%       thetastep = increment in azimuth coordinate theta used as variable 
%          of integration in first surface integral (rad  * TI) - note that
%          it is scaled by 1/T
%       doppbin_frac = presently unused.
%   
% 


gamma = case_params.gammadeg*pi/180;
VR = case_params.VR;
VG = case_params.VG;

delaystep = model_params.acfmodel.delaystep;
maxdelay  = case_params.maxdelay;
acf2matrix  = model_params.acfmodel.acf2matrix;
prn = case_params.prn;
acorrfcn = acf2matrix(:,prn);
acfmatrix  = model_params.acfmodel.acfmatrix;
acfamp = acfmatrix(:,prn);
tauref = model_params.acfmodel.taumatrix;
Ti_input = case_params.Ti;
alt_input = case_params.alt_m;
wdir = case_params.wdirdeg * pi/180;
PDF_params = case_params.PDF_params;
VR = case_params.VR;
VG = case_params.VG;

wavelength = model_params.wavelength;
thetastep  = model_params.thetastep;
%
% This version accounts for the "Doppler Spreading"
%

vs = [ -cos(gamma); 0; -sin(gamma)];
us = [ -cos(gamma); 0; sin(gamma)];


fs = ( dot(VG, vs) - dot(VR, us)) / wavelength;
fc = fs + doppoff;

delay = 1:delaystep:maxdelay;

delay = delay';

ndelay = size(delay, 1);
ref_power = zeros(ndelay,1);

dtheta = thetastep / Ti_input;  % determined experimentally - for
%waveforms
%dtheta = thetastep / (8*Ti_input); % tested for spectrum

theta =  [0:dtheta:2*pi];
ntheta = size(theta,2);

%tic; 

[dp_real, dp_imag, tanx, tany, xo, yo, pdfo, spec_real, fdopp] = ...
   dp_dtheta_dopp_sum_mex(theta, delay, wdir, PDF_params, ...
         Ti_input, alt_input, gamma, 0, 0, fc, VR, VG, ttilde, ftilde, dftilde, wavelength);

%     fprintf(' Time for integration: \n'); toc

pwr = trapz (theta(1,:), dp_real, 2) + i*trapz (theta(1,:), dp_imag, 2);


pwr2 = spec_real(:,1); % sum over theta (only points inside the indicator function)
                       % now done in the mex function

%tic;
%
% Replace with my fft-conv function
%
%yy = delaystep * convnfft(acorrfcn, pwr);

yy = delaystep * conv_fft(acorrfcn, pwr);

%yy = delaystep * conv(acorrfcn, pwr);

%corr_spec = delaystep * convnfft(acorrfcn, pwr2);

corr_spec = delaystep * conv_fft(acorrfcn, pwr2);

%corr_spec = delaystep * conv(acorrfcn, pwr2);

%  fprintf(' Time for convolution (FFT): \n'); toc

%corr_spec_delta = spec_imag;

post_corr_delay = delay(1) + ...
    tauref(1):delaystep:delay(end) + tauref(end);

%
% Generation of the covariance matrix - scattered signal
%

if (nargout > 6)
    nlagspre = size(pwr,1);
    nacf = size(acfamp,1);
    nlags = nlagspre + nacf;
    acfrow = [acfamp', zeros(1,nlagspre)];

    for k=1:nlags
       pwr3 = [zeros(1,k-1), acfrow(1:nlags-k+1)] .* [zeros(1,nacf), pwr'];
    %   size(acfamp)
    %   size(pwr3)
       Rtau(k,:) = delaystep * conv_fft(acfamp, pwr3')';
    end
   Rtau = Rtau(2:end,nacf:end-1);
   
   %
   % Generation of covariance matrix for noise
   %

   nlags_acf_half = (size(model_params.acfmodel.taumatrix,1)-1)/2;

   Rtau0 = zeros(nlags-1,nlags+2*nlags_acf_half);

   for k=1:nlags-1
       Rtau0(k,k:k+2*nlags_acf_half) = acfamp';
   end

   Rtau0 = Rtau0(:,nlags_acf_half+1:nlags+nlags_acf_half-1);
   
end



return






