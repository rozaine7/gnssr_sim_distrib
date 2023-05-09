% 
%
% Modified 4/2022 to prepare for releast to CodeOcean
%
% Modified version to allow measuring CPU time and for the sensitivity
% study with respect to step size in theta.  
%  
% Other than that, this sould be the same simulation at the current
% (10/2015) version of simrun.m
%
% HOWEVER - this has not been checked out - DO NOT RUN THIS VERSION EXCEPT
% FOR THESE BENCHMARKS on speed and step size. 
%
% Set up for final versions of TGARS paper
% hard-code for 4 cases, to get the plots to align properly
%
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

eval(['load ',dirname,'tmid']);  % should eventually move this when the new tmid files are created

eval(['diary ',prefix, 'sim.dia']);



specbin =[6,3,2,2];  % set the bin number to be nearest to specular point - set manually

binno_obs = [1, 3, 4, 5, 7, 8];  % bins to generate histograms

histrange = [1.2 1.2 0.5 0.2];    % range of values for histogram (PDF) plots



                      

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
testcases = [1,2,3,4];

ncases = size(testcases,2);
batchsize = 100; % number of waveforms to simulate for each batch 

CN0min = 38;  % Used for the search of best fit C/N0 based upon noise SNR.
CN0max = 42;
CN0step = 0.2;  

nwf = 100; 
%nic = 1000;
nic = 100;
nicplot = 50;

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
% Flags 
%

stepsize_flag = 0;  % Run the simulator with various step sizes in theta
ensemblesize_flag = 0;   % Run the simuator with various ensemble sizes
plotflag = 1;    % plot results from previoulsy saved results



%
% Begin run of synthetic data generator 
%


if stepsize_flag == 1
%
% Sensitivity of model to the step size. 
%


   fprintf('Simulator sensitivity study...\n')
  
   eval(['load ', prefix,'sim']);
   %
   %  Set up the model parameters (mp) struture
   %
   clear mp

   mp.dftilde = 5;
   mp.dttilde = 5e-4;
   mp.ftilde_max = 600;
   mp.ttilde_max = 0.03;
   
   Bmin = 0;
   

   
   eval(['load ', prefix,'CN0']);

   dtheta_range = 2*pi * [1/4, 1/5, 1/10, 1/20, 1/50, 1/100]*0.001  % number are in fractions of a circle
   
   
   simtime_th = zeros(size(dtheta_range,2), size(testcases,2));
   comptime_th = zeros(size(dtheta_range,2), size(testcases,2));
   
   simcpu_th = zeros(size(dtheta_range,2), size(testcases,2));
   compcpu_th = zeros(size(dtheta_range,2), size(testcases,2));
   
   %
   % Group waveform data ranges - these will be usd for comparison
   %  estimator has already been run. 
   %
   % In this section, the model is run for the 4 cases selected from the data
   % 
   
   binno_syn2 = interp1(delaymean_syn,[1:length(delaymean_syn)], ...
        delaymean_obs(binno_obs),'nearest');
    
   corrspec_bin1 = zeros(size(ftilde,2), size(testcases,2),  size(dtheta_range,2));
   corrspec_bin2 = zeros(size(ftilde,2), size(testcases,2),  size(dtheta_range,2));
   corrspec_bin3 = zeros(size(ftilde,2), size(testcases,2),  size(dtheta_range,2));
   corrspec_bin4 = zeros(size(ftilde,2), size(testcases,2),  size(dtheta_range,2));
   corrspec_bin5 = zeros(size(ftilde,2), size(testcases,2),  size(dtheta_range,2));
   corrspec_bin6 = zeros(size(ftilde,2), size(testcases,2),  size(dtheta_range,2));
   
   ntheta = zeros(size(dtheta_range, 2),1);
           
   for ktheta = 1:size(dtheta_range,2)
       
      for kplot = testcases;


          cp.gammadeg = elhist(kplot); 
          cp.alt_m = althist(kplot);
          cp.VR =  VRscat(kplot,:);
          cp.VG = VGscat(kplot,:); 
          cp.prn = prnhist(kplot);
          cp.roll=0; %rollint(j)*pi/180;
          cp.pitch=0; %pitchint(j)*pi/180;
          cp.yaw=0; %yawint(j)*pi/180-azm(j,ind); for now ...
       
          fprintf('Simulating case %5i: SV%02d at el: %2.2f deg: MSS= %8.6f \n', ...
             kplot, cp.prn,cp.gammadeg, xest(kplot,1) );
         
       %
       % The gold-code ACF is used so this needs to be re-run at each
       % example.  All of the case_params elements are dummies (though once
       % this works, I could re-substitue those given above ...
       %
           [mp, cpdummy] = modelset( mp, 8, 0, 10000, pi/2, 1e-3, 0, ...
              [0 0 0], [0 0 0], 10000, prnhist(kplot), 'CA', dtheta_range(ktheta));  % set up the simuato model parameters
       
           cp.PDF_params = [sqrt(xest(kplot,1)) sqrt(xest(kplot,1)) 0 0  0  0 0 0];

           fprintf('Step size: Delta Theta %8.6e \n', mp.thetastep );
            
           [pcdarray, y2, ftilde, corrspec, fdsurf, Rtau,dummy, Rtau0, ...
            comptime_th(ktheta, kplot), simtime_th(ktheta, kplot),  ...
            compcpu_th(ktheta, kplot), simcpu_th(ktheta, kplot), ntheta(ktheta)] = ...
            wfsim(mp, cp, CN0(kplot), 0, 0, 1, Bmin, 1.0);  % generates synthetic waveforms
    
     
       %
       % Comparison of spectral model
       %
           
           corrspec_bin1(:,kplot, ktheta) = corrspec(binno_syn2(1),:);
           corrspec_bin2(:,kplot, ktheta) = corrspec(binno_syn2(2),:);
           corrspec_bin3(:,kplot, ktheta) = corrspec(binno_syn2(3),:);
           corrspec_bin4(:,kplot, ktheta) = corrspec(binno_syn2(4),:);
           corrspec_bin5(:,kplot, ktheta) = corrspec(binno_syn2(5),:);
           corrspec_bin6(:,kplot, ktheta) = corrspec(binno_syn2(6),:);


   end
   
   end
   
   save benchmarks 
   
end


%
% Timing of ensemble sim vs. computational sim
%

if ensemblesize_flag == 1
     nwf = 100
     ensemble_size = [10, 100, 1000, 1e4];
     
     dtheta_ens = 0.0003;
     
       fprintf('Ensemble size sensitivity study...\n')
  
   eval(['load ', prefix,'sim']);
   %
   %  Set up the model parameters (mp) struture
   %
   clear mp

   mp.dftilde = 5;
   mp.dttilde = 5e-4;
   mp.ftilde_max = 600;
   mp.ttilde_max = 0.03;
   
   Bmin = 0;
   

   
   eval(['load ', prefix,'CN0']);
     
      
   simtime_en = zeros(size(ensemble_size,2), size(testcases,2));
   comptime_en = zeros(size(ensemble_size,2), size(testcases,2));
   
   simcpu_en = zeros(size(ensemble_size,2), size(testcases,2));
   compcpu_en = zeros(size(ensemble_size,2), size(testcases,2));
   
   
   

       for kens = 1:size(ensemble_size,2)
           
              for kplot = testcases;


          cp.gammadeg = elhist(kplot); 
          cp.alt_m = althist(kplot);
          cp.VR =  VRscat(kplot,:);
          cp.VG = VGscat(kplot,:); 
          cp.prn = prnhist(kplot);
          cp.roll=0; %rollint(j)*pi/180;
          cp.pitch=0; %pitchint(j)*pi/180;
          cp.yaw=0; %yawint(j)*pi/180-azm(j,ind); for now ...
       
          fprintf('Simulating case %5i: SV%02d at el: %2.2f deg: MSS= %8.6f \n', ...
             kplot, cp.prn,cp.gammadeg, xest(kplot,1) );
         

       [mp, cpdummy] = modelset( mp, 8, 0, 10000, pi/2, 1e-3, 0, ...
              [0 0 0], [0 0 0], 10000, prnhist(kplot), 'CA', dtheta_ens);  % set up the simuato model parameters
       
       cp.PDF_params = [sqrt(xest(kplot,1)) sqrt(xest(kplot,1)) 0 0  0  0 0 0];
          
           fprintf('Ensemble size: %8.6e \n', ensemble_size(kens) );
            
           [pcdarray, y2, ftilde, corrspec, fdsurf, Rtau,dummy, Rtau0, ...
            comptime_en(kens, kplot), simtime_en(kens, kplot),  ...
            compcpu_en(kens, kplot), simcpu_en(kens, kplot), nthetaens] = ...
            wfsim(mp, cp, CN0(kplot), nwf, ensemble_size(kens), 1, Bmin, 1.0);  % generates synthetic waveforms
         
            fprintf('Ensemble size check:  \n' );
            size(y2)
            fprintf('step size check:  \n' );
            ntheta
            
       end
       end
       
       save benchmarks

end


if plotflag == 1
    %
    % Plotting 
    %
   
    load benchmarks 
    
       plotcase = 1;  % set which of the simulated cases is to be plotted. 
   plottheta = [2, 3, 4, 5, 6];  % which step size results to plot
   

  %
  % Plot of the spectra for comparison
  %
  %set(0,'DefaultAxesFontName', 'Times New Roman')
  set(0,'DefaultAxesFontSize', 12)

  % Change default text fonts.
   %set(0,'DefaultTextFontname', 'Times New Roman')
  set(0,'DefaultTextFontSize', 10)

  
  figure(1)
  %
  % Bin 2
  %
  subplot(2,2,1)
  plot(ftilde', squeeze(corrspec_bin2(:,1,plottheta(1))), 'k', ...
       ftilde', squeeze(corrspec_bin2(:,1,plottheta(2))), 'k--', ...
        ftilde', squeeze(corrspec_bin2(:,1,plottheta(3))), 'b', ...
        ftilde', squeeze(corrspec_bin2(:,1,plottheta(4))), 'b--', ...
        ftilde', squeeze(corrspec_bin2(:,1,plottheta(5))), 'r')
  xlabel('f-tilde(Hz)')
  ylabel('S(f-tilde)')
  title(['Lag = ', ...
            num2str(round(10*pcdarray(binno_syn2(2))/293)/10),' chips'])
  hl1 = legend(['n = ',num2str(ntheta(plottheta(1)))], ...
      ['n = ',num2str(ntheta(plottheta(2)))], ...
      ['n = ',num2str(ntheta(plottheta(3)))], ...
      ['n = ',num2str(ntheta(plottheta(4)))], ...
    ['n = ',num2str(ntheta(plottheta(5)))], 'Location','NorthWest');


  set(hl1, 'FontSize',  9)

  axis([-600 600 0 3e-3])
  
  %
  % Bin 4
  %
  subplot(2,2,2)
  
  plot(ftilde', squeeze(corrspec_bin4(:,1,plottheta(1))), 'k', ...
       ftilde', squeeze(corrspec_bin4(:,1,plottheta(2))), 'k--', ...
        ftilde', squeeze(corrspec_bin4(:,1,plottheta(3))), 'b', ...
        ftilde', squeeze(corrspec_bin4(:,1,plottheta(4))), 'b--', ...
        ftilde', squeeze(corrspec_bin4(:,1,plottheta(5))), 'r')
  xlabel('f-tilde(Hz)')
  ylabel('S(f-tilde)')
  title(['Lag = ', ...
            num2str(round(10*pcdarray(binno_syn2(4))/293)/10),' chips'])
%  legend(['n = ',num2str(ntheta(plottheta(1)))], ...
%      ['n = ',num2str(ntheta(plottheta(2)))], ...
%      ['n = ',num2str(ntheta(plottheta(3)))], ...
%      ['n = ',num2str(ntheta(plottheta(4)))], ...
%    ['n = ',num2str(ntheta(plottheta(5)))]);
  
    axis([-600 600 0 2e-3])
    
  %
  % Bin 5
  %
  subplot(2,2,3)
  
  plot(ftilde', squeeze(corrspec_bin5(:,1,plottheta(1))), 'k', ...
       ftilde', squeeze(corrspec_bin5(:,1,plottheta(2))), 'k--', ...
        ftilde', squeeze(corrspec_bin5(:,1,plottheta(3))), 'b', ...
               ftilde', squeeze(corrspec_bin5(:,1,plottheta(4))), 'b--', ...
        ftilde', squeeze(corrspec_bin5(:,1,plottheta(5))), 'r')
  xlabel('f-tilde(Hz)')
  ylabel('S(f-tilde)')
  title(['Lag = ', ...
            num2str(round(10*pcdarray(binno_syn2(5))/293)/10),' chips'])
  %legend(['n = ',num2str(ntheta(plottheta(1)))], ...
  %    ['n = ',num2str(ntheta(plottheta(2)))], ...
  %    ['n = ',num2str(ntheta(plottheta(3)))], ...
  %    ['n = ',num2str(ntheta(plottheta(4)))], ...
  %  ['n = ',num2str(ntheta(plottheta(5)))]);

      axis([-600 600 0 5e-4])
  
  % 
  % Bin 6
  %
  subplot(2,2,4)
   
  plot(ftilde', squeeze(corrspec_bin6(:,1,plottheta(1))), 'k', ...
       ftilde', squeeze(corrspec_bin6(:,1,plottheta(2))), 'k--', ...
        ftilde', squeeze(corrspec_bin6(:,1,plottheta(3))), 'b', ...
               ftilde', squeeze(corrspec_bin6(:,1,plottheta(4))), 'b--' , ...
        ftilde', squeeze(corrspec_bin6(:,1,plottheta(5))), 'r')
  xlabel('f-tilde(Hz)')
  ylabel('S(f-tilde)')
  title(['Lag = ', ...
            num2str(round(10*pcdarray(binno_syn2(6))/293)/10),' chips'])
 % legend(['n = ',num2str(ntheta(plottheta(1)))], ...
 %     ['n = ',num2str(ntheta(plottheta(2)))], ...
 %     ['n = ',num2str(ntheta(plottheta(3)))], ...
 %     ['n = ',num2str(ntheta(plottheta(4)))], ...
 %   ['n = ',num2str(ntheta(plottheta(5)))]);
 
       axis([-600 600 0 3e-4])
   
  set(gcf, 'PaperPosition', [0 0 7 5]);  
  set(gcf, 'PaperSize', [7 5])
  
  saveas(gcf, [latexfiledir, 'spectral_error'], 'pdf')
  
  %
  % Plot of runtime vs. step size (probalby not used) and 
  %
  figure(2)
  
  subplot(2,1,1)
  
  plot(ntheta, compcpu_th(:,1), '-ok', ntheta, comptime_th(:,1), '-xk', ...
      ntheta, compcpu_th(:,2:end), '-ok', ntheta, comptime_th(:,2:end), '-xk')
  xlabel('N')
  ylabel('Time (sec)')
  legend('CPU time', 'Stopwatch time', 'Location', 'NorthWest')
  
  timerate_th = comptime_th./(ntheta*ones(1,4))
  cpurate_th = compcpu_th./(ntheta*ones(1,4))
  
  subplot(2,1,2)
      
  loglog(ensemble_size', compcpu_en(:,1), '-.ok', ensemble_size', comptime_en(:,1), '-.xk', ...
      ensemble_size', simcpu_en(:,1), '-ok', ensemble_size', simtime_en(:,1), '-xk', ...
      ensemble_size', compcpu_en(:,2:end), '--ok', ensemble_size', comptime_en(:,2:end), '--xk', ...
      ensemble_size', simcpu_en(:,2:end), '-ok', ensemble_size', simtime_en(:,2:end), '-xk')
       
  xlabel('Ensemble Size')
  ylabel('Time (sec)')
  legend('Model Generation: CPU time', 'Model Generation: Stopwatch time', ...
      'Simulation: CPU time', 'Simulation: Stopwatch time', 'Location', 'NorthWest')
  
  timerate_th = simtime_en./(ensemble_size'*ones(1,4))
  cpurate_th = simcpu_en./(ensemble_size'*ones(1,4))
  
  set(gcf, 'PaperPosition', [0 0 8 5]);  
  set(gcf, 'PaperSize', [8 5])

  saveas(gcf, [latexfiledir, 'timing'], 'pdf')
  
end




diary off



