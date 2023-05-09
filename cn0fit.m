% 
% Set up for final versions of TGARS paper
% hard-code for 4 cases, to get the plots to align properly
%
% Minor fixes JLG 08/13
% First distribution for CYGNSS project 10/12
% 
%  Added new consistency checks - changes the order of bin-bin covariance
%  matrix calcualtion
%
%  Presently configured to read an example run from an airborne retrieval
%  "level2" file. 
%
%  Clean up for IGARSS paper.  JLG 5/12
% 
% Deleted the consistency tests from modeltest_data.m - will load the data
% from an estimation run and saved as ***level2.mat [for now hardcoded]
%
% Script to test model simulations.  First comparison - DMR data from PALS
%  2008 underflight
% First written for Venture proposal
%
clear

%
% load the saved files - give the directory in "dirname" and the prefix in
%  "prefix".  File must have name <dirname><prefix>level2.mat
%
dirname = './';
prefix = '2Mar09_'


fprintf(['Loading ', prefix,'level2.mat ... \n'])

eval(['load ',dirname, prefix,'level2']);

eval(['diary ',prefix, 'sim.dia']);

%
% Diretory of latex files for paper - where all generated plots will go
%
latexfiledir = '/Users/jgarriso/jgarriso/pubs/papers/2012_tgars_simulator/draft1/'; 

specbin =[6,3,2,2];  % set the bin number to be nearest to specular point - set manually

binno_obs = [1, 3, 4, 5, 7, 8];  % bins to generate histograms

histrange = [1.2 1.2 0.5 0.2];    % range of values for histogram (PDF) plots



%
% These flags control the script so that you don't need to re-run the
% simulator in order to re-generate plots.
%
% Make sure that these are not saved into a file that will be re-loaded!
%
% master_run_flag presently unused.
%

master_run_flag  = 0;  % 1=process data from scratch
                       % 2=process data from scratch and plot
                       % 0 = plot only from saved data
                       % -1 = do not process data or plot
                       
noise_fit = 1;         % 1= estimate C/N0 from noise floor bin.                       
                       
master_sim_flag = 0;    % 1 = run simumations 
                       % 0 = don't run simulation
                       
master_plot_flag = 0;    % 1 = plot results
                       % 0 = do not plot results. 
                       
model_check = 0;       % 1= perform internal consistency check - will verify
                       %    synthetic data agrees with numerical model in
                       %    both time and frequency domain.
                       % 0 = skip the internal consistency check.
                       
synthetic_level1B = 0;  % 1= generate a synthetic level1B file with the simulated
                        %    waveforms to be re-run throught the estimator.
                        %    file name:  <prefix>sim_level1B.mat
                        % Not impelemnted as of 1/31/2015
                        
                        
                      

%
% These options set how the simulator will run. 
%
% testcases = indicies in the level2 batches - the simulator will be
% independently run at the conditions of each one of these batches. 
%
% nic waveforms will be generated - each one incoherent averaging of nwf
% waveforms.
% Integraton time set in cp.Ti
%
%testcases = [1,2,3,4];
testcases = 1;
ncases = size(testcases,2);
batchsize = 100; % number of waveforms to simulate for each batch 
CN0 = power(10,[3.9, 3.9, 4.0, 4.0]); % right now, the reflectivity of H2O must be taken into 
                     % account here. 
%CN0 = power(10, 3.9)

CN0min = 38;  % Used for the search of best fit C/N0 based upon noise SNR.
CN0max = 42;
CN0step = 0.2;

nwf = 100; 
%nic = 1000;
nic = 1000;

fignum =1;


%
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
% C/N0 best fit
%

%if (noise_fit == 1)

   fprintf('Running simulator ...\n')
   
   %
   %  Set up the model parameters (mp) struture
   %
   clear mp

   mp.dftilde = 5;
   mp.dttilde = 5e-4;
   mp.ftilde_max = 600;
   mp.ttilde_max = 0.03;
   
   Bmin = 100;

   %
   % Group waveform data ranges - these will be usd for comparison
   %  estimator has already been run. 
   %
   % In this section, the model is run for the 4 cases selected from the data
   % 
   
   CN0dBtest = CN0min:CN0step:CN0max; 
   ndB = size(CN0dBtest,2);
   
   PDFobs = zeros(ncases, 10);
   PDFsyn = zeros(ncases, ndB, 10);
   PDFargobs = zeros(ncases, 10);
   PDFargsyn = zeros(ncases, ndB, 10);

   for kdB = 1:ndB  % loop over carrier to noise ratio.
      CN0test = power(10, CN0dBtest(kdB)/10);
      for kplot = testcases
       
          cp.gammadeg = elhist(kplot); 
          cp.alt_m = althist(kplot);
          cp.VR =  VRscat(kplot,:);
          cp.VG = VGscat(kplot,:); 
          cp.prn = prnhist(kplot);
          cp.roll=0; %rollint(j)*pi/180;
          cp.pitch=0; %pitchint(j)*pi/180;
          cp.yaw=0; %yawint(j)*pi/180-azm(j,ind); for now ...
       
          fprintf('Simulating case %5i: SV%02d at el: %2.2f deg: MSS= %8.6f CN0 =%10.4f dB\n', ...
              kplot, cp.prn,cp.gammadeg, xest(kplot,1), CN0dBtest(kdB) );
       %
       % The gold-code ACF is used so this needs to be re-run at each
       % example.  All of the case_params elements are dummies (though once
       % this works, I could re-substitute those given above ...
       %
         [mp, cpdummy] = modelset( mp, 8, 0, 10000, pi/2, 1e-3, 0, ...
            [0 0 0], [0 0 0], 10000, prnhist(kplot), 'CA');  % set up the simuato model parameters
         
          cp.PDF_params = [sqrt(xest(kplot,1)) sqrt(xest(kplot,1)) 0 0  0  0 0 0];

          [pcdarray, y2, ftilde, corrspec, fdsurf, Rtau,dummy, Rtau0] = ...
             wfsim(mp, cp, CN0test, nwf, nic, 1, Bmin, 1.0);  % generates synthetic waveforms
     
          wfsyn = y2;
          delsyn = pcdarray/293;    % in chips
          
          [wfnobs, nfobs, delaymean_obs, totalpwr_obs]  = wfnorm( ...
             del_inp((kplot-1)*nwf_use(1)+1:kplot*nwf_use(kplot),:)-xest(kplot,2), ...
             wf_inp((kplot-1)*nwf_use(1)+1:kplot*nwf_use(kplot),:) ); 
   
          [wfnsyn, nfsyn, delaymean_syn, totalpwr_syn]  = wfnorm( ...
             squeeze(delsyn(1,:)), squeeze(wfsyn(:,:)) );
          varobs = var(wfnobs(:,1));
          varsyn = var(wfnsyn(:,1));
          vardiff = varobs-varsyn;
          fprintf('  Noise power:  Obs = %10.7f  Syn = %10.7f  Diff =%10.7f \n', ...
              varobs, varsyn, vardiff );

          
          varobsmat(kplot,kdB) = varobs;
          varsynmat(kplot,kdB) = varsyn;
          vardiffmat(kplot,kdB) = vardiff;
          
          [Nobs, histbin_obs] = hist(wfnobs(:,1));
          [Nsyn, histbin_syn] = hist(wfnsyn(:,1));
          totalarea_obs = trapz(histbin_obs, Nobs);
          totalarea_syn = trapz(histbin_syn, Nsyn);
          PDFobs(kplot,:) = Nobs/(totalarea_obs);
          PDFsyn(kplot, kdB, :) = Nsyn/(totalarea_syn);
          PDFargobs(kplot,:) = histbin_obs;
          PDFargsyn(kplot, kdB, :) = histbin_syn;
  
    
      end
      
    
      %%%%%
   
   end
   