function [y, ysignal, yn, fdata, finterp, specinterp, H, W] = datasim(f, spec, ts, n, m, CN0, powervariation, Rbin,Rtau0)
%
% Mod 10/12 JLG - change order of bin-bin correlation to implimnet this
% first before the time-correlation filtering. 
%
% simulate random data, of n points, with sample period of ts sec.  
% Data will be normally distributed, with a power spectrum of spec(f),
% where f is in Hz.  spec has dimensions of nfreq * nbins 
%
% m sets of the n points will be incoherently averaged.  If m=0 then only
% one set of the n points will be returned (no averageing).
%
% CN0 = carrier to noise ratio in W-Hz.
%
% Rbin = cross-correlation between bins - must be PSD -
  
% 
% convert the spectrum to a discrete frequency response function.
%
tdata = ts * n;  % length of data
fs = 1/tdata;   % frequency resolution
fdata = 1/ts;    % Nyquist rate.

%P = fdata;  % Power of incoming, discrete time white noise signal with
%A = sqrt(P); % PSD of 1 (pass through anti-aliasing filter with BW of fdata). 

A = powervariation;  % incoming signal has a power of 1 W. 

if CN0 == 0
  PN = 0;
else
  PN = 1/(CN0 * ts);  % Power in noise given CN0 and defining C==1
end
AN = sqrt(PN/2);    % simulate a complex signal with PN 

nbins = size(spec,2);

finterp = -fdata/2 + [0:(n-1)]' * fs;
specinterp = interp1(f, spec, finterp);

specinterp(logical(specinterp<0))=0;  % interpolation will sometimes give a negative value

H = sqrt(specinterp);

%
% Bin-bin correlation matrix Rbin - must be PSD
%
if nargin > 7
   W = chol(Rbin);
   
 %  [V0, D0] = eig(Rtau0);  % very crude attempt to stop numerical problems wiht chol()
 %  V0 = real(V0);
 %  D0 = real(D0);
   
  %   W0 = D0*sqrt(V0);
  
  W0 = chol(Rtau0);
%   normscale = 1./wf_array';
else
   W = 1;
   W0 = 1;
%   normscale = ones(1,nbins);
end

%
% Loop for each block of n waveforms.
%
if m==0
    % ysignal = A * randcorr(n, nbins, H) * sqrt(fdata);  % randcorr scaling by fdata
    %  ysignal = ysignal * W; 
    %H = -1; % test 
     ysignal = A * randcorr(n, nbins, H, W) * sqrt(fdata);
     yn = AN * (randn(n, nbins) + j * randn(n, nbins));
     y = ysignal + yn;
else
  y = zeros(m,nbins);
  for k=1:m
  
      ycohsig = A *  randcorr(n, nbins, H, W) * sqrt(fdata);  % randcorr scaling by fdata
     % ycohsig = ycohsig .* (ones(n,1) * normscale);  Should be taken care
     % of with the unity power spectrum
 %     ycohsig = ycohsig * W; 
 %    yn = AN * (randn(n, nbins) + j * randn(n, nbins)); % no bin-bin for thermal noise

      yn = AN * randcorr(n, nbins, -1, W0);  % Add bin-bin correlation for thermal noise. 
      ycoh = ycohsig + yn;
      y(k,:) = mean(ycoh .* conj(ycoh), 1);
      
      fprintf(' Generated block number   %5i \n', k)
  end
end

return



function y = randcorr(n, m, H, W)
% 
% Generate m independent columns of correlated time series.  The time
% serisl along the n dimension will be white noise passed through a
% filter with a frequency response function of H. This function is all
% discrete time, the dimensions of H must be set properly before calling.
%
%
% Generate the white noise input sequence and FT it
%
   
   xwhite_indep = (1/sqrt(2)) * randn(n,m) + i * (1/sqrt(2)) * randn(n,m);
   xwhite = xwhite_indep * W;

%
% Filter it using the H(f) from above
%
   if (H ~=-1)
      Xf = fft(xwhite, [], 1);
      Xshift = fftshift(Xf,1);

      Y = H .* Xshift;

      Yshift = ifftshift(Y,1);
      y = ifft(Yshift,[],1);
   else
      y = xwhite;
   end
   
return