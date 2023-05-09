%
% Demonstration of single run of simulator - all input parameters set
% manually.
% 
clear

%
% These flags control the script so that you don't need to re-run the
% simulator in order to re-generate plots.
%
% Make sure that these are not saved into a file that will be re-loaded!
%
% master_run_flag presently unused.
%


%
% These options set how the simulator will run. 

%
% nic waveforms will be generated - each one incoherent averaging of nwf
% waveforms.
%
nwf = 100; % 100 is typical for airbonre

nic = 1000;

% Key data extracted from file:
% File has been georeferenced, due to time limit, will work with this and
%   re-interpolate back to average delay bin location
% 
%   ta0 = time tag (sec) for midpoint of each batch retreival 
%   nwf_use = number of waveforms in each batch - don't know how to handle
%       a varyable number of waveforms - so will assume it is the same for
%       the entire run and only use nwf_use(1).
%   del_inp, wf_inp = delay bins and waveform, in rows of nwf_use each and 
%       nbins columns.  Delay in units of chips. 
%   prnhist = satellite PRN at midpoint of each batch
%   elhist  = satellite elevantion (deg) at midpoint of each batch
%   althist = satellite altitude (m) at midpoint of each batch
%   azhist = satellite azimuth (deg) at midpoint of each batch
%   VGscat = satellite velocity in scattering plane components (m/s)
%   VRscat = receiver veloicty in scattering plane components (m/s) 
%
% Outputs from simulator:
%   wfsyn = waveform (Y^2) in 3-D array, (sample, delay, case)
%   delsyn = corresponding lag in chips.
%


%
%  Set up the model parameters (mp) struture
%

mp.dftilde = 5;  %  step in f-tilde, spectrum independnet variable (Hz)
mp.dttilde = 5e-4; % step i t-tilde, autocorrelation independent variable (sec)
mp.ftilde_max = 600;  % spectrum will be generated for a range of f-tilde from +/- mp.ftilde_max
mp.ttilde_max = 0.03; % autocorrelation will be generated for a range of t-tilde from +/- mp.ttilde_max

%
% Set up the case parameters (cp) structure
%

cp.wdirdeg = 0; % wind direction (deg) - relative to scattering plane
cp.Ti=0.001; % coherent integration time (sec)
cp.nDoppbins = 0; % number of  Doppler bins (one sided)
cp.maxdelay = 2000; % maximum delay to generate waveform (meters)

 
       
%
% airborne case (PALS flight)
%
cp.alt_m = 3.4153e+03;  % altitude (meters)
cp.VR =  [ 12.0212  172.5443   -0.2845];  % receiver velocity (m/s) in 
cp.VG = [ 761.4    2526.0   -855.1] 
cp.gammadeg = 66.4309;
       
%
       % Orbiting case
       %
   %    cp.alt_m = 600e3;
   %    cp.VR =  [ 7500  0  0];
   %    cp.VG = [ 761.4    2526.0   -855.1] 
       
       cp.prn = 15;
       cp.roll=0; %rollint(j)*pi/180;
       cp.pitch=0; %pitchint(j)*pi/180;
       cp.yaw=0; %yawint(j)*pi/180-azm(j,ind); for now ...
       
       CN0 = power(10,4.1);
       Bmin = 0;  % minimum bandwidth of spectrum 
       
   
       %
       % The gold-code ACF is used so this needs to be re-run at each
       % example.  All of the case_params elements are dummies (though once
       % this works, I could re-substitue those given above ...
       %
       [mp, cpdummy] = modelset( mp, 8, 0, 10000, pi/2, 1e-3, 0, ...
           [0 0 0], [0 0 0], 10000, cp.prn, 'CA');  % set up the simuato model parameters
       
       cp.PDF_params = [0.1195 0.1195 0 0  0  0 0 0];

       [pcdarray, y2, ftilde, corrspec, fdsurf, Rtau,dummy, Rtau0] = ...
           wfsim(mp, cp, CN0, nwf, nic, 1, Bmin, 1.0);  % generates synthetic waveforms
    
 
 


