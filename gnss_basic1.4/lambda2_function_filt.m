function [tau, Rf, fplot, frfplot] = ...
      l2filt(BW, dtau, indexmax, n2pts, prn)
%
% 10/22/03 - Flips the Rf about the zero-delay point 
% BW is the single-sided bandwidth
%
% Uses the ifft to compute the autocorrelation waveform for the C/A 
% code signal, after passing it through a filter defined by the 
% frequency response function frf(f), where f ranges from 0 to fmax.
%
% frf should be a column vector
%
% Input units f= Hz, output units, tau=sec
%
% prn input will allow use of the true ACF for the C/A gold codes
%  set prn <= 0 to use ideal Maximal length code.
%
Tc = 9.767e-7; %sec
fc = 1.023e6; % Hz
butterorder = 4; % order of butterworth filter model

nsamples = 2^n2pts;
tau = 0:dtau: (nsamples * dtau);
if (prn <= 0)
   [R, R2] = lambda2_function(tau);
else
   R = ca_acf( prn, tau);
end
df = 1/(max(tau));
fmax = 1/dtau;

%f = [0:df:fmax] * fc; replace in favor of butterworth filter
%
% compute the unfiltered and filtered PSD's
%
psd0 = real(fft(R));
f = [-fmax/2:df:fmax/2]* fc;
fshift = ifftshift(f);
%frf = 1.0 * (f < BW) + 1.0 * ((max(f)-f) < BW);
frf = 1./(1 + (fshift./BW).^(2*butterorder));
psd_filt = psd0 .* frf;

%
% Use ifft to compute the auto-correlations
%
rtau = 2.0 * real(ifft(psd_filt)); % 2 accounts for fact that +/- freq
                                 % included
%Rf = [flipud(rtau(1:indexmax)'); rtau(1:indexmax)'];
Rf = [rtau(nsamples-indexmax+1:nsamples)'; rtau(1:indexmax)'];
tau = [-flipud(tau(1:indexmax)'); tau(1:indexmax)'];

%frfplot = [flipud(frf(1:nsamples/2)'); frf(1:nsamples/2)'];
frfplot = ifftshift(frf);
%fplot = [-flipud(f(1:nsamples/2)'); f(1:nsamples/2)'];
fplot = f;

return