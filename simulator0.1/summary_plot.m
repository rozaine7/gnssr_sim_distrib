%
% Script to print the results of MC retrieval simulations
%
meanplot = 0;  %set to 1 to use mean of the estimates over all baches.


mss = 0.0055:0.002:0.015;
nmss = size(mss,2)
   
tauplot = -1.1:0.01:6;  % in chips - tau axis for plotting best fit
datastructsim.Winv = eye(size(tauplot,2), size(tauplot,2));
datastructsim.wfm  = zeros(size(tauplot));
datastructsim.taum = tauplot;
 
for k=1:nmss

   CN0dB = 45;
 
   subdirname = ['./C',num2str(CN0dB),'MSS',num2str(k),'B'];
%   if  ~exist(subdirname)
%      eval(['mkdir ', subdirname])
%      fprintf(['Created subdirectory ', subdirname,'\n']);   
%   end
   
   eval(['load ',subdirname,'/Level2MSS']);
   
   fprintf(['Loading simulation from ', subdirname,'\n'])
   
   if (meanplot == 1)  % generate a single plot of best fit for whole batch
      mssret = mean(pdfest_array);
      stdmss = std(pdfest_array);
      scalef = mean(scalef_array);
      tau0 = mean(tau_array);
 
      xstate(1) = mssret;
      xstate(2) = tau0;
      xstate(3) = scalef; 
   
      wfmodel  = fwf(xstate, datastructsim, cp, mp);
   
      plot(tauaxis/mp.acfmodel.chipsize, wf, '.b', ...
          datastruct.taum, datastruct.wfm+datastruct.nf, '.r', ...
          tauplot, -wfmodel+datastruct.nf, '-k')
    
      axis([-2 6 0 1.5])
      xlabel('Lag (chips)')
      ylabel('Reflected Waveform (1W input)')
      title([date,' CN0 =',num2str(CN0dB),' MSS =', num2str(mss(k))])
      
      pause
    
      eval(['print wfplotC', num2str(CN0dB),'MSS', ...
        num2str(mss(k)),'.pdf -dpdf'])
   else
      
      nbatches = size(dpinterval,1);
      for kbatch = 1:nbatches
          xstate(1) = pdfest_array(kbatch);
          xstate(2) = tau_array(kbatch);
          xstate(3) = scalef_array(kbatch);
   
          wfmodel  = fwf(xstate, datastructsim, cp, mp);
   
          plot(tauaxis(dpinterval(kbatch,1):dpinterval(kbatch,2),:)/mp.acfmodel.chipsize, ...
            wf(dpinterval(kbatch,1):dpinterval(kbatch,2),:), '.b', ...
            tauplot, -wfmodel+datastruct.nf, '-k')
    
          axis([-2 6 0 1.5])
          xlabel('Lag (chips)')
          ylabel('Reflected Waveform (1W input)')
          title([date,' CN0 =',num2str(CN0dB),' MSS =', num2str(mss(k)), ...
              ' Batch No. ', num2str(kbatch)])
    
          eval(['print wfplotC', num2str(CN0dB),'MSS', ...
            num2str(mss(k)),'b',num2str(kbatch),'.pdf -dpdf'])
           
          pause
      end
   end
      mssarray(k,:) = pdfest_array';
      postcovarray(k,:) = pdf_std';
   
end