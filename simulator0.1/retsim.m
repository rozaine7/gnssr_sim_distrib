%
% Simulated retrievals
%
%eval(['load ', dirname,'/level1B']);
%eval(['load ', dirname, '/runfile ']);

%
% set up parameters for simulated data
%

%taumask = 1:50:1200;  % 1/2 chip spacing
%
% Basic setup parameters
%
load runfile
%cpmp_base

CN0dB = 45;
nwf=100;  % number of waveforms generated and incoherently averaged per sample
nic = 5000; % number of samples
Ti=1e-3;  % coherent integration time

reflectivity = 0.63;  % reflectivity of sea water

%mss = 0.002:0.002:0.02 % Range of MSS to simulate. 
mss = 0.0055:0.002:0.008;
nmss = size(mss,2)

for k=1:nmss

   fprintf('Generating synthetic data for MSS = %6.3f ...\n', mss(k))
   cp.PDF_params = [sqrt(mss(k)) sqrt(mss(k)) 0 0 0 0 0 0];
   [pcd, wf, ftilde, corrspec] = wfsim(mp, cp, power(10,CN0dB/10)*reflectivity, nwf, nic);
   nwf = size(wf,1);
   tauaxis = ones(nwf,1) * pcd;
   prn = 1 * ones(nwf,1);
   twf = [0:nic-1]*Ti*nwf; 
   twf=twf';
   EL = cp.gammadeg * ones(nwf,1);
   altac = cp.alt_m * ones(nwf,1);
   vr_vec = ones(nwf,1) * cp.VR;
   vs_vec = ones(nwf,1) * cp.VG;
   
   subdirname = ['./C',num2str(CN0dB),'MSS',num2str(k),'B'];
   if  ~exist(subdirname)
      eval(['mkdir ', subdirname])
      fprintf(['Created subdirectory ', subdirname,'\n']);   
   end
   
   eval(['!cp runfile.mat ',subdirname]); % always copy a new runfile.mat even if the subdir still exists.
   
   fprintf(['Saving simulated data in ', subdirname,'\n'])
   eval(['save ', subdirname,'/level1B nwf tauaxis wf prn twf EL altac vr_vec vs_vec '])
   
   fprintf('Running estimator \n')
   run_pdf_est(subdirname);
   
end

