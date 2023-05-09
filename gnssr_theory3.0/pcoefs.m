function model_params = pcoefs( model_params)
%
% 1/04 made into a function to use structures.
% Generate the coefficients for a series approximation of the no-Doppler
% waveform
%

model_params.acfmodel = acfgen('CA', 0);

model_params.wavelength = 0.19;
model_params.thetastep  = 1e-5;
model_params.doppbin_frac = 1;


case_params.wdirdeg = 0;
case_params.Ti = 1e-3;
case_params.VR = [0 0 0];
case_params.VG = [0 0 0];
case_params.prn = 1;
case_params.nDoppbins =0;
case_params.alt_m = 10000;
case_params.maxdelay = 20000;

ELref = [20:20:60,70:5:90]';
sqrtMSS = [0.04:0.02:0.15, 0.18:0.06: 0.4];
a1 = zeros(size(ELref,1),size(sqrtMSS,2));
a2 = zeros(size(ELref,1),size(sqrtMSS,2));
a3 = zeros(size(ELref,1),size(sqrtMSS,2));
a4 = zeros(size(ELref,1),size(sqrtMSS,2));
a5 = zeros(size(ELref,1),size(sqrtMSS,2));
a6 = zeros(size(ELref,1),size(sqrtMSS,2));
a7 = zeros(size(ELref,1),size(sqrtMSS,2));
a8 = zeros(size(ELref,1),size(sqrtMSS,2));
a9 = zeros(size(ELref,1),size(sqrtMSS,2));

for q = 1:size(sqrtMSS,2)
  case_params.PDF_params = [sqrtMSS(q), sqrtMSS(q), 0, 0, 0, 0, 0, 0];
  fprintf(' Computing waveforms for sigma = %15.5f\n', sqrtMSS(q));
  for k=1:size(ELref,1)
    case_params.gammadeg = ELref(k);
    [wf, cd, df, pa, fs] = d2map( case_params, model_params);
    cdarray = (cd(2)-cd(1))*[0:1:size(pa,1)-1]';
    [P, S] = polyfit( cdarray/10000, log(10000*pa),8);
    a1(k,q) = P(1);
    a2(k,q) = P(2);
    a3(k,q) = P(3);
    a4(k,q) = P(4);
    a5(k,q) = P(5);
    a6(k,q) = P(6);
    a7(k,q) = P(7);
    a8(k,q) = P(8);
    a9(k,q) = P(9);
  end
end

model_params.series.a1 = a1;
model_params.series.a2 = a2;
model_params.series.a3 = a3;
model_params.series.a4 = a4;
model_params.series.a5 = a5;
model_params.series.a6 = a6;
model_params.series.a7 = a7;
model_params.series.a8 = a8;
model_params.series.a9 = a9;
model_params.series.ELref = ELref;
model_params.series.sqrtMSS = sqrtMSS;

ds = diff(sqrtMSS);
de = diff(ELref);
da1ds = diff(a1,1,2)./(ones(size(a1,1),1)*ds);
da2ds = diff(a2,1,2)./(ones(size(a2,1),1)*ds);
da3ds = diff(a3,1,2)./(ones(size(a3,1),1)*ds);
da4ds = diff(a4,1,2)./(ones(size(a4,1),1)*ds);
da5ds = diff(a5,1,2)./(ones(size(a5,1),1)*ds);
da6ds = diff(a6,1,2)./(ones(size(a6,1),1)*ds);
da7ds = diff(a7,1,2)./(ones(size(a7,1),1)*ds);
da8ds = diff(a8,1,2)./(ones(size(a8,1),1)*ds);
da9ds = diff(a9,1,2)./(ones(size(a9,1),1)*ds);

model_params.series.da1ds = da1ds;
model_params.series.da2ds = da2ds;
model_params.series.da3ds = da3ds;
model_params.series.da4ds = da4ds;
model_params.series.da5ds = da5ds;
model_params.series.da6ds = da6ds;
model_params.series.da7ds = da7ds;
model_params.series.da8ds = da8ds;
model_params.series.da9ds = da9ds;
model_params.series.mids = (sqrtMSS(1:end-1)+sqrtMSS(2:end))/2;

return