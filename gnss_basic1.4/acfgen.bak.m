function acfstruct = acfgen( code_type, prn, FEBW )
  acfstruct = struct('acf2matrix', 0, 'acfmatrix', 0, 'taumatrix', 0, ...
		     'chipsize',293,'delaystep',1);
%
% 1/04 modified to use structures to eliminate global variables.
% complete set of ACF's saved in acfstruct.  All of the component
% functions calls were moved internal to this file.  
%
% Generates an array of code auto-correlation functions (ACFs) for
% the code type specified in code_type, for the satellite defined by 
% prn.
%
% Split from d2map.m code to avoid re-generation each time it is called
% 2/03
%
% Right now, call with a single prn - to generate the ACF's for 1 SV
% set prn=0  and it will generate a matrix, and make it global 
%  

if nargin == 1
  prn = 0;
end

if nargin == 2
  FEBW = 10.23e6;
end

fprintf('Generating the Autocorrelation functions: \n');

if  code_type == 'CA'
  fprintf ('  Infinitely long code \n')
  chipsize = 293;
  delaystep = 3;
  tauref = [-150:150]' * delaystep;  
  [acorr_fcn_volt, acorr_fcn] = lambda2_function (tauref / chipsize);
  if (prn == 0)
     acfmatrix = acorr_fcn_volt * ones(1,32);
  end
elseif  code_type == 'GC' % Gold code ACF's used (no filtering)
  fprintf ('  C/A Gold code, infinite bandwidth \n')
  chipsize = 293;
  delaystep = 3;
  tauref = [-1500:1500]' * delaystep;  
  if (prn == 0)
     acfmatrix = zeros(size(tauref,1),32);
     for kprn = 1:32
	acfmatrix(:,kprn) = ca_acf (kprn, tauref / chipsize);
     end
  else
     acfmatrix = ca_acf (prn, tauref / chipsize);
  end       
elseif  code_type == 'CF' % front-end filtering with ML code
  chipsize = 293;
%  delaystep = 0.3;
  delaystep = 0.0005;
  taumax = 5;
%  taumax = 200;
  fprintf ('  Infinitely long code with %10.1f MHz bandwidth \n', FEBW)
  indexmax = ceil(taumax*(chipsize/delaystep));
  powerof2 = ceil(log2(2*indexmax));
  [tauref_chips, volt_acorr_fcn, fq, PSD] = ...
      lambda2_function_filt (FEBW*1e6, delaystep/chipsize, ...
        indexmax, powerof2+2, 0);
  tauref = tauref_chips * chipsize;
  if (prn == 0)
     acfmatrix = volt_acorr_fcn * ones(1,32);
  end  
elseif  code_type == 'GF' % front-end filtering with C/A sidelobes
  chipsize = 293;
  delaystep = 0.03;
  taumax = 10;
  indexmax = ceil(taumax*(chipsize/delaystep));
  powerof2 = ceil(log2(2*indexmax));
  fprintf ('  C/A Gold code with %10.1f MHz bandwidth \n', FEBW)
  if (prn == 0)
     acfmatrix = [];
     for kprn = 1:32
         [tauref_chips, acorr_fcn_volt, fq, PSD] = ...
           lambda2_function_filt (FEBW*1e6, delaystep/chipsize, ...
           indexmax, powerof2+2, kprn);
	 acfmatrix = [acfmatrix, acorr_fcn_volt];
     end
  else
      [tauref_chips, acorr_fcn_volt, fq, PSD] = ...
         lambda2_function_filt (FEBW*1e6, delaystep/chipsize, ...
         indexmax, powerof2+2, prn);
      acfmatrix = acorr_fcn_volt;
  end       
  tauref = tauref_chips * chipsize;
elseif code_type == 'PY'
  fprintf('  P(Y) Code \n')
  chipsize = 29.3;
  delaystep = 0.3;
  tauref = [-150:150]' * delaystep;
  [lambda, acorr_fcn] = lambda2_function (tauref / chipsize);
elseif  code_type == 'BC'
  fprintf ('  BOC Modulation: \n')
  fs = input('    Subcarrier frequency (fs, in 1.023 MHz) ?\n');
  fc = input('    Chipping frequency (fc, in 1.023 MHz) ? \n');
  fprintf('   Generating ACF for a BOC( %3.1f, %3.1f ) signal.\n', fs, fc);
  chipsize = 293;
  delaystep = 3;
  tauref = [-150:150]' * delaystep;  
  [acorr_fcn_volt, acorr_fcn] = boc_acf (tauref / chipsize, fs ,fc);
  acfmatrix = acorr_fcn_volt * ones(1,32);
else
  fprintf('Error --- wrong code type ...\n')
  return
end

acfstruct.acfmatrix = acfmatrix;
acfstruct.acf2matrix = acfmatrix.^2;
acfstruct.taumatrix =tauref;
acfstruct.chipsize = chipsize;
acfstruct.delaystep = delaystep;
%acfstruct.dacf2dtau = ...
%    diff(acfstruct.acf2matrix,1,1)./(diff(tauref) * ones(1,32));

%
% Subfunctions
%

function [lambda, lambda2] = lambda2_function(tau_chips);
%
%  computes the lambda^2 function in units of half-code chips
%
npts = size(tau_chips,1);
lambda = max( ones(npts,1)-abs(tau_chips), zeros(npts,1));
lambda2 = lambda.^2;

function [tau, Rf, fplot, frfplot] = ...
      lambda2_function_filt(BW, dtau, indexmax, n2pts, prn)
%
% 10/22/03 - Flips the Rf about the zero-delay point 
% BW is the single-sided bandwidth
%
% Uses the ifft to compute the autocorrelation waveform for the C/A 
% code signal, after passing it through a filter defined by the 
% frequency response function frf(f), where f ranges from 0 to fmax.
%
% Input units f= Hz, output units, tau=sec
%
% prn input will allow use of the true ACF for the C/A gold codes
%  set prn <= 0 to use ideal Maximal length code.
%
Tc = 9.767e-7; %sec
fc = 1.023e6; % Hz
nsamples = 2^n2pts;
tau = 0:dtau: (nsamples * dtau);
if (prn <= 0)
   [R, R2] = lambda2_function(tau);
else
   R = ca_acf( prn, tau);
end
df = 1/(max(tau));
fmax = 1/dtau;
f = [0:df:fmax] * fc;
%
% compute the unfiltered and filtered PSD's
%
psd0 = real(fft(R));
frf = 1.0 * (f < BW) + 1.0 * ((max(f)-f) < BW);
psd_filt = psd0 .* frf;
%
% Use ifft to compute the auto-correlations
%
rtau = real(ifft(psd_filt)); % 2 accounts for fact that +/- freq
                                   % included
Rf = [rtau(nsamples-indexmax+1:nsamples)'; rtau(1:indexmax)'];
tau = [-flipud(tau(1:indexmax)'); tau(1:indexmax)'];

frfplot = [flipud(frf(1:nsamples/2)'); frf(1:nsamples/2)'];
fplot = [-flipud(f(1:nsamples/2)'); f(1:nsamples/2)'];

function R = ca_acf( prn, tau )
%
%  Generates the ACF of the C/A "Gold codes".  
%

[pnmatrix, G1, G2] = pngen(0);

maxlag = ceil(max(tau));

[Rtable, lags] = xcorr( pnmatrix(:,prn), pnmatrix(:,prn), ...
                            maxlag, 'unbiased');

R = interp1( lags, Rtable, tau);

function [pnmatrix, G1, G2] = pngen(prn);
%
% function to generate the PN codes with amplitudes from -1 to 1
%
% prn is a matrix of satellite numbers (PRN's actually)
% pnmatrix is the output array of dimension 1023 X size(prn)
%
% G1 is the maximal length shift prn code from shift register G1.
% G2 is the full state of the 10 stage shift register G2 (1023 X 10)
%
% Generate the G1 codes.
%
G1 = zeros(1023,1);
G2 = zeros(1023,10);
G1_register = [1 1 1 1 1 1 1 1 1 1];
G2_register = [1 1 1 1 1 1 1 1 1 1];
for G1_counter = 1:1023
  G1(G1_counter,1) = G1_register(10);
  feedback1 = xor(G1_register(3), G1_register(10));
  G1_register = [feedback1, G1_register(1:9)];
end
%
% Generate the G2 codes
%
for G2_counter = 1:1023
  G2(G2_counter,:) = G2_register;
  feedback2 = mod(sum(G2_register([2,3,6,8,9,10])),2);
  G2_register = [feedback2, G2_register(1:9)];
end
%
% Generate each of the 32 PRNS
%

taps1 = [2, 3, 4, 5, 1, 2, 1, 2, 3, 2, 3, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, ...
	 6, 1, 4, 5, 6, 7, 8, 1, 2, 3, 4];
taps2 = [6, 7, 8, 9, 9, 10, 8, 9, 10, 3, 4, 6, 7, 8, 9, 10, 4, 5, 6, 7, ...
	 8, 9, 3, 6, 7, 8, 9, 10, 6, 7, 8, 9];

pnmatrix = 2*mod( (G1*ones(1,32) + G2(:,taps1) + G2(:,taps2)),2)-1;


function [R, R2] = boc_acf( tau, fs, fc)
% 
% ACF for BOC(fs, fc) modulated signals
%  tau is in units of TC   (tau =  tau_sec / Tc_sec) to stay consistent
%  with the other ACF functions. 
% fs, fc are the multiple of 1.023 MHz
%
%  See Betz, 2001 Navigation paper.

n = 2 * fs/fc
ell = -n:1:n;   % index of peaks (2n-1 total)
peak_ts = ell';  % Location of peaks, in units of Ts
peak_tc = peak_ts/n;  % location of peaks in units of Tc.
RY_peak = ((-1).^ell ) .* ( n - abs(ell))/n; % value of peaks

% Values of the ACF at each peak.
peak_tc = [min(tau); min(peak_tc) - 1/n; peak_tc; max(peak_tc) + 1/n; max(tau)];
RY_peak = [0; 0; RY_peak'; 0; 0];

% Interpolated to the lags, tau

R = interp1(peak_tc, RY_peak, tau);
R2 = R.^2;

return






