function [tau, rtau] = filter_waveform(fmax, frf, Tc)
%
% Uses the ifft to compute the autocorrelation waveform for the C/A 
% code signal, after passing it through a filter defined by the 
% frequency response function frf(f), where f ranges from 0 to fmax.
%
% frf should be a column vector
%
% Input units f= Hz, output units, tau=sec
%
%Tc = 9.767e-7; %sec
nsamples = size(frf,1);
df = fmax / (nsamples-1);
f=[0:nsamples-1]' *df;
%
% compute the unfiltered and filtered PSD's
%
psd0 = sinc( Tc * f).^2;
psd_filt = psd0 .* frf;
%
% Use ifft to compute the auto-correlations
%
rtau = real(ifft( psd_filt));
dtau = 1/fmax;
tau = [0:nsamples-1]'*dtau;
return