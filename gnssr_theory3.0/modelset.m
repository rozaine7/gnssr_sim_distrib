function [model_params, case_params] = modelsetup( model_params, ...
    ws_cm, wdir, alt, gamma, Ti, nDoppbins, VR, VG, maxdelay, prn, code_type, thetastep, BW, fc)
%
%  Mod 4/16 by JLG - to accomadate SoOps
%   BW = bandwidth of the SoOp (assuming BPSK for now)
%   fc = carrier frequency
%  Test script to set all of the parameters in the model_params structure
%  Mod 9/08 by JLG to have a function-line input option, in addition to the
%  interactive interface. 
% 
%  [model_params, case_params] = modelset(model_params, windspeed, winddir,
%                            altitude, elevation, coherent_int_time, nDopp,
%                            Vreceiver, VGPS, maxdelay, PRN, code_type)
%
model_params.c = 2.99792458e8;  %m/s

if (nargin <= 1)
     ws_cm = input('Wind Speed (for Cox and Munk Model) ');
   case_params.wdirdeg = input('Wind Direction (deg) ');
   case_params.alt_m = input('Receiver altitude (m) ');
   case_params.gammadeg = input('Elevation (deg) ');
   case_params.Ti = input('Predetection Integration Time (sec) ');
   case_params.nDoppbins = input('Number of doppler bins to create ');
     VR(1) = input('Receiver velocity, i component (m/s)');
     VR(2) = input('Receiver velocity, j component (m/s)');
     VR(3) = input('Receiver velocity, k component (m/s)');
     VG(1) = input('GPS velocity, i component (m/s)');
     VG(2) = input('GPS velocity, j component (m/s)');
     VG(3) = input('GPS velocity, k component (m/s)');
   case_params.maxdelay = input('Max delay to compute range bins (m) ');
   case_params.prn = input('PRN Code (only for Gold Code option ');
   code_type = input('ACF Model CA=ideal, GC=Gold Code, SO = Signal of opportunity )');
   if code_type == 'SO'
       BW = input(' Signal of opportunity transmitted bandwidth ');
       fc = input(' Carrier frequency (MHz)');
   else
       BW = -999;
   end
else
   case_params.wdirdeg = wdir;
   case_params.alt_m = alt;
   case_params.gammadeg = gamma;
   case_params.Ti = Ti;
   case_params.nDoppbins = nDoppbins;
   case_params.maxdelay = maxdelay;
   case_params.prn = prn;
   if code_type ~= 'SO'
       BW = -999;
   end
end
   case_params.VR = VR;
   case_params.VG = VG;
   if code_type ~= 'SO'
      coxmunk_params = coxmunk_pdf (ws_cm);
      PDF_params(1) = sqrt(coxmunk_params(1));
      PDF_params(2) = sqrt(coxmunk_params(2));
      PDF_params(3:7) = coxmunk_params(3:7);
      PDF_params(8) = 0;
      case_params.PDF_params = PDF_params;
   
      model_params.wavelength = 0.19;
   else
       model_params.c
       fc
      model_params.wavelength = model_params.c / (fc * 1e6);  % meters, fc in MHz. 
      kstar = 2*pi/(3 * model_params.wavelength); 
      [sigma2_u, sigma2_c, sigma2omni] = var_from_spec(ws_cm, 0.84, kstar)
      PDF_params(1) = sqrt(sigma2_u);
      PDF_params(2) = sqrt(sigma2_c);
      PDF_params(3:8) = 0;

      case_params.PDF_params = PDF_params;
   end
   fullmodel.PDF = 'gram_charlier';
 

   model_params.acfmodel = acfgen( code_type, 0, 10e6, BW);

   model_params.fullmodel = fullmodel;
   model_params.doppbin_frac = 1;
   
   if code_type ~= 'SO'
      
   else 
 
   end
   if nargin > 12
      model_params.thetastep  = thetastep;
   else
      model_params.thetastep  = 1e-5;
   end
return
