function [wfi, dYdtaui, dYdsigmai ] =  pswf( tau, case_params, model_params)
%
% 10/04 - included gradient option.  Will return dY/dM and dY/tau as
%  column vectors, if they are asked for in the return
%  argument. M=sqrtMSS^2 = the MSS
% 1/04 - made to use the new structures
% Extract approximation to waveform shape from powerseries at lags tau
%
% Input tau is in meters
%
  
acf2matrix = model_params.acfmodel.acf2matrix;
dacf2dtau   = model_params.acfmodel.dacf2dtau;
taumatrix  = model_params.acfmodel.taumatrix;
chipsize   = model_params.acfmodel.chipsize;
delaystep  = model_params.acfmodel.delaystep;

a1 = model_params.series.a1;
a2 = model_params.series.a2;
a3 = model_params.series.a3;
a4 = model_params.series.a4;
a5 = model_params.series.a5;
a6 = model_params.series.a6;
a7 = model_params.series.a7;
a8 = model_params.series.a8;
a9 = model_params.series.a9;
ELref = model_params.series.ELref;
sqrtMSS = model_params.series.sqrtMSS;

alt = case_params.alt_m;
sigma0 = case_params.PDF_params(1);
satel = case_params.gammadeg;

%
% Interpolating the series coefficients to the correct Elevation and
% MSS
%
a1i = interp2( sqrtMSS, ELref, a1, sigma0, satel);
a2i = interp2( sqrtMSS, ELref, a2, sigma0, satel);
a3i = interp2( sqrtMSS, ELref, a3, sigma0, satel);
a4i = interp2( sqrtMSS, ELref, a4, sigma0, satel);
a5i = interp2( sqrtMSS, ELref, a5, sigma0, satel);
a6i = interp2( sqrtMSS, ELref, a6, sigma0, satel);
a7i = interp2( sqrtMSS, ELref, a7, sigma0, satel);
a8i = interp2( sqrtMSS, ELref, a8, sigma0, satel);
a9i = interp2( sqrtMSS, ELref, a9, sigma0, satel);

P = [a1i, a2i, a3i, a4i, a5i, a6i, a7i, a8i, a9i];


dt = taumatrix(2)-taumatrix(1);  % takes tau spacing from that saved with acf
nlambda = size(taumatrix,1);
taurange = max(max(tau)) - min(min(tau));
taup = 0:dt:taurange+600;  % found that it sometimes is necessary to
                                % run the power out this far to avoid
                                % numerical problem
acorrfcn = acf2matrix(:,1);
dacorrfcn = dacf2dtau(:,1);
taubar = taup / alt;
pabar = exp(polyval( P, taubar));
pa = pabar / alt;
wf = conv(acorrfcn, pa);
wfnorm = wf / (sum(wf,2) * (dt/chipsize));
%wfnorm = wf / (dt/chipsize);
tauwf = -dt*floor(nlambda/2) + dt * [1:size(wf,2)];

padding = -2000:10:-dt*nlambda-10;

tauwf = [padding, tauwf];
wfnorm = [zeros(size(padding)), wfnorm];
wfi = interp1(tauwf',wfnorm',tau);

%
% gradient calculation
%

if nargout > 1
  %
  % Computed at the same steps as pabar (tauwf).
  %
  
  da1ds = model_params.series.da1ds;
  da2ds = model_params.series.da2ds;
  da3ds = model_params.series.da3ds;
  da4ds = model_params.series.da4ds;
  da5ds = model_params.series.da5ds;
  da6ds = model_params.series.da6ds;
  da7ds = model_params.series.da7ds;
  da8ds = model_params.series.da8ds;
  da9ds = model_params.series.da9ds;
  mids = model_params.series.mids;
  
  da1i = interp2( mids, ELref, da1ds, sigma0, satel);
  da2i = interp2( mids, ELref, da2ds, sigma0, satel);
  da3i = interp2( mids, ELref, da3ds, sigma0, satel);
  da4i = interp2( mids, ELref, da4ds, sigma0, satel);
  da5i = interp2( mids, ELref, da5ds, sigma0, satel);
  da6i = interp2( mids, ELref, da6ds, sigma0, satel);
  da7i = interp2( mids, ELref, da7ds, sigma0, satel);
  da8i = interp2( mids, ELref, da8ds, sigma0, satel);
  da9i = interp2( mids, ELref, da9ds, sigma0, satel);

  dads = [da1i, da2i, da3i, da4i, da5i, da6i, da7i, da8i, da9i];
  dpbardsigma      = pabar .* polyval( dads, taubar);
  
  dpdsigma = dpbardsigma/ alt;
    
  dYdtau = conv(dacorrfcn, pa) / (sum(wf,2) * dt/chipsize);
  dYdsigma = conv(acorrfcn, dpdsigma) / (sum(wf,2) * dt/chipsize);

  dYdtau = [zeros(size(padding)), dYdtau];
  dYdsigma = [zeros(size(padding)), dYdsigma];

  dYdtaui = interp1(tauwf(2:end)', dYdtau', tau);
  dYdsigmai = interp1(tauwf', dYdsigma', tau);
  
end

return 
