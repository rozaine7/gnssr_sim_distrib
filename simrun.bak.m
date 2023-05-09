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

eval(['load ',dirname,'tmid']);  % should eventually move this when the new tmid files are created

eval(['diary ',prefix, 'sim.dia']);



specbin =[6,3,2,2];  % set the bin number to be nearest to specular point - set manually

binno_obs = [1, 3, 4, 5, 7, 8];  % bins to generate histograms

histrange = [1.2 1.2 0.5 0.2];    % range of values for histogram (PDF) plots

dtheta = 0.0013;

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
                       
noise_fit = 0;         % 1= estimate C/N0 from noise floor bin. 
                       % -1 = use pre-defined C/N0 values hardcoded at that
                       % point.
                       % 0 or anthing else - use the C/N0 values 
                       
master_sim_flag = 0;    % 1 = run simumations 
                       % 0 = don't run simulation
                       
master_plot_flag = 1;    % 1 = plot results
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
testcases = [1,2,3,4];

ncases = size(testcases,2);
batchsize = 100; % number of waveforms to simulate for each batch 

CN0min = 38;  % Used for the search of best fit C/N0 based upon noise SNR.
CN0max = 42;
CN0step = 0.2;  

nwf = 100; 
%nic = 1000;
nic = 1000;
nicplot = 50;

dtheta = 0.001;

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

if (noise_fit == 1)

   fprintf('CN0 estimator ...\n')

   %
   %  Set up the model parameters (mp) struture
   %
   clear mp

   mp.dftilde = 5;
   mp.dttilde = 5e-4;
   mp.ftilde_max = 600;
   mp.ttilde_max = 0.03;
   
   Bmin = 0;

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
            [0 0 0], [0 0 0], 10000, prnhist(kplot), 'CA', dtheta);  % set up the simuato model parameters
         
          cp.PDF_params = [sqrt(xest(kplot,1)) sqrt(xest(kplot,1)) 0 0  0  0 0 0];

          [pcdarray, y2, ftilde, corrspec, fdsurf, Rtau,dummy, Rtau0] = ...
             wfsim(mp, cp, CN0test, nwf, nic, 1, Bmin, 1.0);  % generates synthetic waveforms
     
          wfsyn = y2;
          delsyn = pcdarray/293;    % in chips
          
          [wfnobs, nfobs, delaymean_obs, totalpwr_obs]  = wfnorm( ...
             del_inp((kplot-1)*nwf_use(kplot)+1:kplot*nwf_use(kplot),:)-xest(kplot,2), ...
             wf_inp((kplot-1)*nwf_use(kplot)+1:kplot*nwf_use(kplot),:) ); 
   
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
   
   
   end
   
   
   [mindiff, indexmindiff] = min(abs(vardiffmat'));
   CN0 = power(10, CN0dBtest(indexmindiff)/10)
   
   fprintf('Saving best-fit CN0  ...\n')
   eval(['save ', prefix,'CN0'])

elseif (noise_fit == -1)
    
%
% hardcoded C/N0 values
%
   fprintf('Using pre-defined CN0  ...\n')
   CN0 = power(10,[3.9, 3.9, 4.0, 4.0]); % right now, the reflectivity of H2O must be taken into 
                     % account here. 
   eval(['save ', prefix,'CN0']);
end

%%%%%
%%%%%
%
% Model Consistency Check
%
if (model_check == 1)

   fprintf('Performing consistency check ...\n')
   %
   %  Set up the model parameters (mp) struture
   %
   clear mp

   mp.dftilde = 5;   % right now these are hard-coded - should the be moved ?
   mp.dttilde = 5e-4;
   mp.ftilde_max = 600;
   mp.ttilde_max = 0.03;
          
   binno_check = [150, 200, 250, 450];  % bins for comparison in the consistency check
   
   eval(['load ', prefix,'CN0']);

   %
   % In this section, the model is run for the 4 cases selected from the data
   % 

   for kplot = testcases;

       cp.gammadeg = elhist(kplot);  % Elevation (deg)
       cp.alt_m = althist(kplot);    % Altitude (md)
       cp.VR =  VRscat(kplot,:);    % Receiver velocity - in scattering plane coordinates
       cp.VG = VGscat(kplot,:);     % GPS velocity - in scattering plane coordinates
       cp.prn = prnhist(kplot);     % GPS PRN 
       cp.roll=0; %rollint(j)*pi/180;  % aircraft attitude - if used. Not all cases saved it
       cp.pitch=0; %pitchint(j)*pi/180;
       cp.yaw=0; %yawint(j)*pi/180-azm(j,ind); for now ...

       
       fprintf('Simultor consitency check for case %5i: SV%02d at el: %2.2f deg \n', ...
           kplot, cp.prn,cp.gammadeg );
       %
       % The gold-code ACF is used so this needs to be re-run at each
       % example.  All of the case_params elements are dummies (though once
       % this works, I could re-substitue those given above ...
       %
       [mp, cpdummy] = modelset( mp, 8, 0, 10000, pi/2, 1e-3, 0, ...
           [0 0 0], [0 0 0], 10000, prnhist(kplot), 'CA', dtheta);
       
       cp.PDF_params = [sqrt(xest(kplot,1)) sqrt(xest(kplot,1)) 0 0  0  0 0 0];

       % 
       % Generate synthetic complex waveforms. 
       %
  
       [pcdarray, dummy1, dummy2, dummy3, dummy4, Rtau_raw] = ...
              wfsim(mp, cp, 0, 2000, 0, -1, 0, 1.0);
       [pcdarray, y_bbcorr1, ftilde, corrspec, fdsurf, Rtau] = ...
              wfsim(mp, cp, 0, 2000, 0, 1, 0, 1.0);
       [pcdarray, y_bbcorr0, ftilde, corrspec, fdsurf, Rtau, wf_array] = ...
              wfsim(mp, cp, 0, 2000, 0, 0, 0, 1.0);
       
       wf_bbcorr1 = mean(y_bbcorr1 .* conj(y_bbcorr1),1);
       wf_bbcorr0 = mean(y_bbcorr0 .* conj(y_bbcorr0),1);
       
       wf_spec = trapz(ftilde, corrspec, 2);
       
       figure(fignum)
       
       plot(pcdarray, wf_array, pcdarray, wf_bbcorr1, 'o', ...
            pcdarray, wf_bbcorr0, 'x', pcdarray, diag(Rtau_raw), ...
            pcdarray, wf_spec, pcdarray(binno_check), wf_spec(binno_check), '*')
       legend('Integrated', 'Simulated, R(tau,tau)', ...
           'Simulated, R(tau,tau)=I', 'diag(R(tau,tau))', '\int S(f) df', ...
           'Delay bins for spectral test')
       axis([-500 1500 -0.1 1.1])
       xlabel('Delay (m)')
       ylabel('Y^2(\tau)')
       
       fignum= fignum+1;
       
       %
       % Compare the coherence time for the time-domain and
       % frequency-domain models
       %
       fprintf('Computing time-domain model ... \n')
       
       ttilde = [-mp.ttilde_max:mp.dttilde:mp.ttilde_max];  % sec
       nttilde = size(ttilde,2);

       for k=1:size(ttilde,2)
          [yy(:,k), pcdarray, dummy, pwr, fs, fdsurf ] = ...
                wf_from_gcp_dopp( 0, cp, mp, ttilde(k), 0, mp.dftilde);
          fprintf('   ttilde = %12.5f ms \n', ttilde(k)*1000)
       end

       ftilde_fft = -(1/(2*mp.dttilde)):(1/(2*mp.ttilde_max)):(1/(2*mp.dttilde));
       corrspec_fft = fftshift(fft(ifftshift(yy),[],2))*mp.dttilde;

       figure(fignum)
       for ksubplot = 1:4
           [Sf1, fwelch] = pwelch(y_bbcorr1(:,binno_check(ksubplot)),[],[],[],1/(cp.Ti));
           [Sf0, fwelch] = pwelch(y_bbcorr0(:,binno_check(ksubplot)),[],[],[],1/(cp.Ti));

           subplot(2,2,ksubplot)
           plot(ftilde, corrspec(binno_check(ksubplot),:), ...
               ftilde_fft, real(corrspec_fft(binno_check(ksubplot),:)), ...
               fwelch-max(fwelch)/2, fftshift(Sf1), 'x', ...
               fwelch-max(fwelch)/2, fftshift(Sf0), 'o')
           xlabel('f-tilde (Hz)')
           ylabel('S(f)')
           legend('S(f) model', 'FFT(R(\tau))', 'Numerical: R(\tau, \tau)', ...
               'Numerical: R(\tau,\tau) = I')
           axis([-500 500 -min(Sf1) max(Sf1)*1.1]) 
       end
      % pause
    
   end
   
    eval(['save ', prefix,'sim']);

end

%
% Begin run of synthetic data generator 
%

if (master_sim_flag == 1)

   fprintf('Running simulator ...\n')
   
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

   %
   % Group waveform data ranges - these will be usd for comparison
   %  estimator has already been run. 
   %
   % In this section, the model is run for the 4 cases selected from the data
   % 

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
           [0 0 0], [0 0 0], 10000, prnhist(kplot), 'CA', dtheta);  % set up the simuato model parameters
       
       cp.PDF_params = [sqrt(xest(kplot,1)) sqrt(xest(kplot,1)) 0 0  0  0 0 0];

       [pcdarray, y2, ftilde, corrspec, fdsurf, Rtau,dummy, Rtau0] = ...
           wfsim(mp, cp, CN0(kplot), nwf, nic, 1, Bmin, 1.0);  % generates synthetic waveforms
    
       wfsyn(:,:, kplot) = y2;
       delsyn(1,:,kplot) = pcdarray/293;    % in chips
    
   end
   
    eval(['save ', prefix,'sim']);

end

if (master_plot_flag ==1 ) 
%
% generate comparison plots for simulated data
% Assumes that the experimental data has been processed and the simulations
% run
%
eval(['load ', prefix,'sim']);

%
% Diretory of latex files for paper - where all generated plots will go
%
%latexfiledir = '/Users/jgarriso/jgarriso/pubs/papers/2012_tgars_simulator/revision/'; 
latexfiledir = '/Users/jgarriso/jgarriso/projects/gnssr/airborne_simulator_paper/airborne_sim_100115/'


subplotid = {'(a)', '(b)', '(c)', '(d)'};

% Note - these are setup-commands, so they only take place on restarting.
%  I need to figure out more how matlab handles this, but for now I'll keep
%  them here. 
% Change default axes fonts.
%set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

% Change default text fonts.
%set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 12)

%
% figure 5 \label{fig:delaycoord} in latex file - demonstration that the
% delay coordinate does not exceed 
%

figure(5)
subplotno = 1;

for kplot = testcases
    subplot(2,2,subplotno)
    hourstart = floor(twf(t_inds(kplot)-floor(nwf_use(1)/2))/3600)
    plot( twf(t_inds(kplot)-floor(nwf_use(1)/2):t_inds(kplot)+floor(nwf_use(1)/2)-1)-hourstart*3600,...
    del_inp((kplot-1)*nwf_use(1)+1:kplot*nwf_use(kplot),specbin(kplot))*mp.acfmodel.chipsize, 'k') 
    title(subplotid(kplot))
    xlabel(['Seconds past ', num2str(hourstart), 'H'])
    ylabel('Spec. Pt. Bin Delay (m)')
    subplotno = subplotno+1;
end

% attempt to set and scale figure to be a direct drop-in for the latex file

set(gcf, 'PaperPosition', [0 0 6 5]);  
set(gcf, 'PaperSize', [6 5])

saveas(gcf, [latexfiledir, 'delaycoord'], 'pdf')

% pause

%
% figure 6 - \label{fig:simvstheory_wf} in latex file - illustration of
% simulated vs. measured waveforms
%
figure(6)
subplotno = 1;

offset = zeros(size(testcases,2), 6);

for kplot = testcases

     [wfnobs, nfobs, delaymean_obs, totalpwr_obs]  = wfnorm( ...
         del_inp((kplot-1)*nwf_use(kplot)+1:kplot*nwf_use(kplot),:)-xest(kplot,2), ...
         wf_inp((kplot-1)*nwf_use(kplot)+1:kplot*nwf_use(kplot),:) ); 
   
     [wfnsyn, nfsyn, delaymean_syn, totalpwr_syn]  = wfnorm( ...
         squeeze(delsyn(1,:,kplot)), squeeze(wfsyn(1:nicplot,:,kplot)) );
     
     binno_syn = interp1(delaymean_syn,[1:length(delaymean_syn)], ...
        delaymean_obs(binno_obs),'nearest');
    
     %
     % Compute best fit from estimator results
     %
     
     cp.gammadeg = elhist(kplot);  % Elevation (deg)
     cp.alt_m = althist(kplot);    % Altitude (md)
     cp.VR =  VRscat(kplot,:);    % Receiver velocity - in scattering plane coordinates
     cp.VG = VGscat(kplot,:);     % GPS velocity - in scattering plane coordinates
     cp.prn = prnhist(kplot);     % GPS PRN 
     cp.roll=0; %rollint(j)*pi/180;  % aircraft attitude - if used. Not all cases saved it
     cp.pitch=0; %pitchint(j)*pi/180;
     cp.yaw=0; %yawint(j)*pi/180-azm(j,ind); for now ...

     [mp, cpdummy] = modelset( mp, 8, 0, 10000, pi/2, 1e-3, 0, ...
           [0 0 0], [0 0 0], 10000, prnhist(kplot), 'CA', dtheta);
       
     cp.PDF_params = [sqrt(xest(kplot,1)) sqrt(xest(kplot,1)) 0 0  0  0 0 0];
  
     [wf_bestfit, chips_bestfit, df, pwr, fs] = d2map( cp, mp );
       
     [wf_bestfitn, nfbestfit, delaymean_bestfitn, totalpwr_bestfit]  = wfnorm( ...
         chips_bestfit'/mp.acfmodel.chipsize, wf_bestfit' );
     
     subplot(4,2,subplotno)
     %plot( del_inp((kplot-1)*nwf_use(1)+1:kplot*nwf_use(kplot),:)-xest(kplot,2), wfnobs, '.k', ...
     %    chips_bestfit/mp.acfmodel.chipsize, xest(kplot,3)*wf_bestfit, ...
     %    [delaymean_obs(binno_obs); delaymean_obs(binno_obs)], ...
     %     [-0.1*ones(size(binno_obs)); 1.6*ones(size(binno_obs))], '--k', ...
     %     delaymean_syn, mean(wfnsyn,1),'xr')
     
     
     plot( del_inp((kplot-1)*nwf_use(1)+1:kplot*nwf_use(kplot),:)-xest(kplot,2), wfnobs, '.k', ...
         delaymean_bestfitn, wf_bestfitn, 'k',...
         [delaymean_obs(binno_obs); delaymean_obs(binno_obs)], ...
          [-0.1*ones(size(binno_obs)); 1.6*ones(size(binno_obs))], '--k')

      
      chips_bestfit_chips = chips_bestfit'/mp.acfmodel.chipsize;
      
      binno_fit = interp1(chips_bestfit_chips,[1:length(chips_bestfit)], ...
        delaymean_obs(binno_obs(2:6)),'nearest');
      
      offset(kplot,1) = 0; 
      offset(kplot,2:6) = xest(kplot,3)*wf_bestfit(binno_fit)'- mean(wfnobs(:,binno_obs(2:6)),1);
      
   %   delaymean_obs(binno_obs)
  
     title(strcat(subplotid(kplot), 'Measured'))     
 %    legend('Measured - 1/2 Chip sampling')
  %   xlabel(['Delay (chips)'])
     axis([-1.7 2.7 -0.1 2.1])
     subplotno = subplotno+1
     
     
     subplot(4,2,subplotno)
     plot( delaymean_syn, wfnsyn,'.k', ...
                  [delaymean_syn(binno_syn); delaymean_syn(binno_syn)], ...
          [-0.1*ones(size(binno_syn)); 1.6*ones(size(binno_syn))], '--k');
      
     title(strcat(subplotid(kplot), 'Simulated'))
   %  xlabel(['Delay (chips)'])
     axis([-1.7 2.7 -0.1 2.1])
     subplotno = subplotno+1

%    xlabel('Lag')
%    ylabel('Waveform (normalized)')
%    title(['C/N0 ', num2str(10*log10(CN0)), 'dB-Hz'])
  clear wfnsyn
  clear delaymean_syn
end

% attempt to set and scale figure to be a direct drop-in for the latex file

set(gcf, 'PaperPosition', [0 0 6 8]);  
set(gcf, 'PaperSize', [6 8])

saveas(gcf, [latexfiledir, 'simvstheory_wf'], 'pdf')

%pause


%
% figure 7 - \label{simvstheory_pdfa-d} in latex file - comparison of
% measured vs. synthetic PDF's, numerically computed as histograms
%
fignum = 7;


for kplot = testcases
    %
    % plot of mean waveform, simulated one, and model - this part is
    % redundnat
    %
    %  subplot(1,2,2)
    % plot( delaymean_syn, mean(wfnsyn,1), delaymean_obs, mean(wfnobs,1), 'x');
    % legend('Synthetic','Measured')
    % axis([-1 6 -0.1 1.2])
    % xlabel('Lag')
    % ylabel('Waveform (normalized)')
    %  title(['C/N0 ', num2str(10*log10(CN0)), 'dB-Hz'])
    
    [wfnobs, nfobs, delaymean_obs, totalpwr_obs]  = wfnorm( ...
         del_inp((kplot-1)*nwf_use(kplot)+1:kplot*nwf_use(kplot),:)-xest(kplot,2), ...
         wf_inp((kplot-1)*nwf_use(kplot)+1:kplot*nwf_use(kplot),:) ); 
   
    [wfnsyn, nfsyn, delaymean_syn, totalpwr_syn]  = wfnorm( ...
         squeeze(delsyn(1,:,kplot)), squeeze(wfsyn(:,:,kplot)) );
  
    binno_syn = interp1(delaymean_syn,[1:length(delaymean_syn)], ...
        delaymean_obs(binno_obs),'nearest');
    
    fprintf([' Noise floor variance (obs): ', num2str(var(wfnobs(:,binno_obs(1)))),'\n'])
    fprintf([' Noise floor variance (obs): ', num2str(var(wfsyn(:,binno_syn(1)))),'\n'])
    fprintf([' Noise floor variance (obs): ', ...
        num2str(var(wfsyn(:,binno_syn(1)))-var(wfnobs(:,binno_obs(1)))),'\n'])
    
    figure(fignum);
    
    for ksubplot=1:6
        subplot(3,2,ksubplot)
        [Nobs, histbin_obs] = hist(wfnobs(:,binno_obs(ksubplot)));
        [Nsyn, histbin_syn] = hist(wfnsyn(:,binno_syn(ksubplot)));
        totalarea_obs = trapz(histbin_obs, Nobs);
        totalarea_syn = trapz(histbin_syn, Nsyn);
       %
       % Showed that mean of wfnsyn is almost identical to the "best fit" 
       %
        mean_obs = mean(wfnobs(:,binno_obs(ksubplot)));
        mean_syn = mean(wfnsyn(:,binno_syn(ksubplot)));
        offsetmean = mean_syn - mean_obs
     %   plot(histbin_obs, Nobs/sum(Nobs), 'x', histbin_syn, Nsyn/sum(Nsyn), 'o');
         plot(histbin_obs, Nobs/(totalarea_obs), '-ok', ...
             histbin_syn, Nsyn/(totalarea_syn), '--b', ...
             histbin_syn-offsetmean, Nsyn/(totalarea_syn), '-r');
        if ksubplot==1
           legend('Observed', 'Synthetic','Synthetic - unbiased')
        end
           
        xlabel('Z')
        ylabel('p_{Y^2}(Z)')
    %    axis([-0.1 1.2 0 max(Nobs/(totalarea_obs))])
        title(['Lag = ', ...
            num2str(round(10*delaymean_obs(binno_obs(ksubplot)))/10),' chips'])
      % round off the delay to first decimal place in title
      %...
      %         'C/N0 ', num2str(10*log10(CN0(kplot))), 'dB-Hz']);    

    end

    fprintf(['Got here - fig number', num2str(fignum),'\n'])      
    set(gcf, 'PaperPosition', [0 0 6 8]);  
    set(gcf, 'PaperSize', [6 8])

    saveas(gcf, [latexfiledir, 'simvstheory_pdf',num2str(kplot)], 'pdf')
    
          
    set(gcf, 'PaperPosition', [0 0 6 8]);  
    set(gcf, 'PaperSize', [6 8])

    %pause

   % fignum = fignum+1;
end


%
% figure 11 - \label{ref:simvstheory_cov} - covariance matrix shown as
% filled countour plots
%

    figure(11)
    subplotno = 1;

fignum = 14

for kplot=testcases
    
    %
    % Selected rows of covariance matrix
    %
    

        figure(fignum)
    
      for ksubplot=1:6
           subplot(3,2,ksubplot)

    %    plot(delaymean_obs, Pobs(binno_obs(ksubplot),:)/Pobs(binno_obs(ksubplot),binno_obs(ksubplot)), 'o', ...
    %       delaymean_syn, Psyn(binno_syn(ksubplot),:)/Psyn(binno_syn(ksubplot),binno_syn(ksubplot)) )
    
          [wfnobs, nfobs, delaymean_obs, totalpwr_obs]  = wfnorm( ...
            del_inp((kplot-1)*nwf_use(kplot)+1:kplot*nwf_use(kplot),:)-xest(kplot,2), ...
            wf_inp((kplot-1)*nwf_use(kplot)+1:kplot*nwf_use(kplot),:) ); 
   
          [wfnsyn, nfsyn, delaymean_syn, totalpwr_syn]  = wfnorm( ...
             squeeze(delsyn(1,:,kplot)), squeeze(wfsyn(:,:,kplot)) );
     
    
           Pobsn = corrcoef(wfnobs);
           Psynn = corrcoef(wfnsyn);

           plot(delaymean_obs, Pobsn(binno_obs(ksubplot),:), 'ok', ...
                delaymean_syn, Psynn(binno_syn(ksubplot),:), '-k' );
           if ksubplot==1
              legend('Observed', 'Synthetic')
           end
           
           ylabel('C(\tau,:)')
           xlabel('Delay (chips)')
           
           axis([-2 5 -0.2 1.2])
     
           title(['\tau = ', ...
               num2str(round(10*delaymean_obs(binno_obs(ksubplot)))/10),' chips'])
           if ksubplot==1
              legend('Observed', 'Synthetic')
           end
          
        
      end
    
    set(gcf, 'PaperPosition', [0 0 6 8]);  
    set(gcf, 'PaperSize', [6 8])

    saveas(gcf, [latexfiledir, 'simvstheory_cov_row',num2str(kplot)], 'pdf')
    
   % pause

   %fignum = fignum+1;
    
end

%%%
%
%  Wind speed and PALS data
%
fignum = fignum + 1;

figure(fignum)

subplot(2,1,1)
[axh, h1, h2] = plotyy(t_pol, wpol, t_pol, tbv)


axis(axh(1), [20 26 0 40])
axis(axh(2), [20 26 120 130])
set(h1,'Color',[0 0 0])
set(h2,'Color',[0 0 1])
set(h2,'LineStyle','--')
set(axh(1), 'YTick', [ 0 10 20 30 40 ])
set(axh(1), 'YColor', [ 0 0 0 ])
set(axh(2), 'YColor', [ 0 0 0 ])

xlabel('Time (Hr)')
ylabel(axh(1),'POLSCAT Wind (m/s)')
ylabel(axh(2),'PALS   T_{b,V} (K)')
legend('POLSCAT', 'PALS')
hold on

bartimes = [tmid/3600  tmid/3600 ]
%plot(t_pol(select), wpol(select), 'o')
plot(bartimes, [0 40], 'k')
text(bartimes(1)+0.1, 5, '(a)');
text(bartimes(2)+0.1, 5, '(b)');
text(bartimes(3)+0.1, 5, '(c)');
text(bartimes(4)+0.1, 15, '(d)');




hold off

subplot(2,1,2)
[axh, h1, h2] = plotyy(twf/3600, prn, twf/3600, EL)
axis(axh(1), [20 26 0 40])
axis(axh(2), [20 26 0 90])
set(h1,'Color',[0 0 0])
set(h2,'Color',[0 0 0])
set(h1,'LineStyle','--')
set(axh(2), 'YTick', [ 0 15 30 45 60 75 90 ])
set(axh(1), 'YTick', [ 0 10 20 30 ])
set(axh(1), 'YColor', [ 0 0 0 ])
set(axh(2), 'YColor', [ 0 0 0 ])

xlabel('Time (Hr)')
ylabel(axh(1),'PRN')
ylabel(axh(2),'Elevation (deg)')
legend('PRN', 'Elevation')

hold on
%plot(t_pol(select), 20, 'o')
plot(bartimes, [0 40], 'k')
text(bartimes(1)+0.1, 5, '(a)');
text(bartimes(2)+0.1, 5, '(b)');
text(bartimes(3)+0.1, 5, '(c)');
text(bartimes(4)+0.1, 5, '(d)');

hold off



fprintf(['Got here - fig number', num2str(fignum),'\n'])      
set(gcf, 'PaperPosition', [0 0 6 8]);  
set(gcf, 'PaperSize', [6 8])

saveas(gcf, [latexfiledir, 'experiment_pdf'], 'pdf')


%%%
eval(['save ', prefix,'sim']);

end

if (synthetic_level1B == 1)
    eval(['load ', prefix,'sim']);
    
    npts = size(wfsyn,1);
    
    subsample_array = delaymean_obs;  % later, we would want to jitter
    subsample_ptr = interp1(delaymean_syn,[1:length(delaymean_syn)], ...
        subsample_array,'nearest'); 
    
    tauaxis = ones(npts,1)* delaymean_syn(subsample_ptr)*mp.acfmodel.chipsize;  
       % in meters - to agree with prepro output.
    wf = squeeze(wfsyn(:,subsample_ptr,1))*4.5e6;  % to make it approximately the same order of mag. 
    
    EL = cp.gammadeg * ones(npts,1);
    AZ = 0 * ones(npts,1);  % fixed for now
    altac = cp.alt_m*ones(npts,1);
    prn = cp.prn *ones(npts,1);
    gpspos = ones(npts,1)* [ 0 0 0];
    recpos = ones(npts,1)* [ 0 0 0];
    gpsvel = ones(npts,1)* cp.VG;
    recvel = ones(npts,1)* cp.VR;
    acatt = ones(npts,1)* [cp.roll, cp.pitch, cp.yaw]; 
    
    eval(['save ',dirname,prefix,'sim_level1B ', ...
      ' tauaxis wf twf EL AZ altac prn latlons gpspos recpos gpsvel recvel acatt']);
end



diary off



