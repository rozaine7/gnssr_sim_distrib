%
% coherent averaging of the simulated data at various integration times.
% to generate plot as requested by TGARS reviewers
%
% Will use the thinned arrays
%

ndelays = size(delsyn_integ001t,2);
shortsum = 0;

%
% TI = 0.002 ms
%
wfsum001 = zeros(100000,ndelays,4);
wfsum002 = zeros(50000,ndelays,4);
wfsum005 = zeros(20000,ndelays,4);
wfsum01 = zeros(10000,ndelays,4);

wfsum001_short = zeros(100000,ndelays,4);
wfsum002_short = zeros(50000,ndelays,4);
wfsum005_short = zeros(20000,ndelays,4);
wfsum01_short = zeros(10000,ndelays,4);

if shortsum == 1
for k=1:100000
  wfsum001_short(k,:,:) = sum(wfsyn_integ0005t(2*k-1:2*k,:,:),1)/2;   
end

for k=1:50000
  wfsum002_short(k,:,:) = sum(wfsyn_integ0005t(4*k-1:4*k,:,:),1)/4;   
end

for k=1:20000
  wfsum005_short(k,:,:) = sum(wfsyn_integ0005t(10*k-9:10*k,:,:),1)/10; 
end


for k=1:10000
  wfsum01_short(k,:,:) = sum(wfsyn_integ0005t(20*k-19:20*k,:,:),1)/20; 
end

end
    

%for k=1:100000
%  wfsum001(k,:,:) = sum(wfsyn_integ001t(2*k-1:2*k,:,:),1)/2;   
%end

for k=1:50000
  wfsum002(k,:,:) = sum(wfsyn_integ001t(2*k-1:2*k,:,:),1)/2;   
end

for k=1:33333
  wfsum003(k,:,:) = sum(wfsyn_integ001t(3*k-2:3*k,:,:),1)/3;   
end

for k=1:20000
  wfsum005(k,:,:) = sum(wfsyn_integ001t(5*k-4:5*k,:,:),1)/5; 
end


for k=1:10000
  wfsum01(k,:,:) = sum(wfsyn_integ001t(10*k-9:10*k,:,:),1)/10; 
end


%pwr005sum = abs(wfsum005).^2;
%pwr002sum = abs(wfsum002).^2;
%pwr001 = abs(wfsyn_integ001).^2;
%pwr002 = abs(wfsyn_integ002).^2;
%pwr005 = abs(wfsyn_integ005).^2;

%snr001 = std(pwr001,1)./mean(pwr001,1);
%snr002 = std(pwr002,1)./mean(pwr002,1);
%snr005 = std(pwr005,1)./mean(pwr005,1);

%snr002sum = std(pwr002sum,1)./mean(pwr002sum,1);
%snr005sum = std(pwr005sum,1)./mean(pwr005sum,1);

plotlables = ['(a)', '(b)', '(c)', '(d)'];

decmate = 1:2:ndelays;

figure(1)

subplotid = {'(a)', '(b)', '(c)', '(d)'};
%set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

% Change default text fonts.
%set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 12)

axisset = [-1.5 2.5 0 0.5;
           -1.5 2.5 0 0.5;
           -1.5 2.5 0 0.6; 
           -1.5 2.5 0 0.9];

for kplot = 1:4
  subplot(2,2,kplot)

  plot(squeeze(delsyn_integ001t(:,:,kplot)), var(squeeze(real(wfsyn_integ001t(:,:,kplot))),1), ':k', ...
    squeeze(delsyn_integ002t(:,:,kplot)), var(squeeze(real(wfsyn_integ002t(:,:,kplot))),1), '--k', ...
    squeeze(delsyn_integ002t(:,:,kplot)), var(squeeze(real(wfsum002(:,:,kplot))),1), 'ok', ...
    squeeze(delsyn_integ003t(:,:,kplot)), var(squeeze(real(wfsyn_integ003t(:,:,kplot))),1), '-k', ...
    squeeze(delsyn_integ003t(:,:,kplot)), var(squeeze(real(wfsum003(:,:,kplot))),1), 'xk')
  %  squeeze(delsyn_integ005t(:,:,kplot)), var(squeeze(real(wfsyn_integ005t(:,:,kplot))),1), '-k', ...
  %  squeeze(delsyn_integ005t(:,:,kplot)), var(squeeze(real(wfsum005(:,:,kplot))),1), 'xk', ...
  %  squeeze(delsyn_integ01t(:,:,kplot)), var(squeeze(real(wfsyn_integ01t(:,:,kplot))),1), ':k', ...
  %  squeeze(delsyn_integ01t(:,:,kplot)), var(squeeze(real(wfsum01(:,:,kplot))),1), '+k')
    
  xlabel('Lag (chips)')
  ylabel('St. Dev. Real(Y)') 
  title(subplotid(kplot))
    axis(axisset(kplot,:))
end

legend('T_I = 1 ms', 'T_I = 2 ms', 'T_I = 1 ms, N_C=2', 'T_I = 3 ms', 'T_I = 1 ms, N_C=3')

figure(2)

for kplot = 1:4
  subplot(2,2,kplot)

  plot(squeeze(delsyn_integ001t(:,:,kplot)), mean(squeeze(abs(wfsyn_integ001t(:,:,kplot)).^2),1), '--k', ...
    squeeze(delsyn_integ005t(:,:,kplot)), mean(squeeze(abs(wfsyn_integ005t(:,:,kplot)).^2),1), '-k', ...
    squeeze(delsyn_integ01t(:,:,kplot)), mean(squeeze(abs(wfsyn_integ01t(:,:,kplot)).^2),1), ':k', ...
    squeeze(delsyn_integ005t(:,decmate,kplot)), mean(squeeze(abs(wfsum005(:,decmate,kplot)).^2),1), 'ok', ...
    squeeze(delsyn_integ01t(:,decmate,kplot)), mean(squeeze(abs(wfsum01(:,decmate,kplot)).^2),1), 'xk')
    
  xlabel('Lag (chips)')
  ylabel('Y^2')
  legend('T_I = 1 ms', 'T_I = 5 ms', 'T_I = 10 ms', '\Sigma^{2} Y', '\Sigma^{5} Y')
  
end

%
% Histograms
%

binno = [150, 250, 350];
binno = 150;

figure(3)

%for kplot = 1:4
%  subplot(2,2,kplot)

%        [N001, histbin001] = hist(pwr001(:,binno,kplot));
%        [N005, histbin005] = hist(pwr005(:,binno,kplot));
%        [N01, histbin01] = hist(pwr01(:,binno,kplot));
 %       [N005sum, histbin005sum] = hist(pwr005sum(:,binno,kplot));
%        [N01sum, histbin01sum] = hist(pwr01sum(:,binno,kplot));
        
%        totalarea001 = trapz(histbin001, N001);
%        totalarea005 = trapz(histbin002, N005);
%        totalarea01 = trapz(histbin005, N01);
%        totalarea_sum005 = trapz(histbin005sum, N005sum);
%        totalarea_sum01 = trapz(histbin01sum, N01sum);
        
%        plot(histbin001, log(N001/(totalarea001)), '--k', ...
%             histbin005, log(N005/(totalarea005)), '-k', ...
%             histbin01, log(N01/(totalarea01)), ':k', ...
%             histbin005sum, log(N005sum/(totalarea_sum005)), 'ok', ...
%             histbin01sum, log(N01sum/(totalarea_sum01)), 'xk')
% ylabel('ln p(Y)')
%  legend('T_I = 1 ms', 'T_I = 2 ms', 'T_I = 5 ms', '\Sigma^{2} Y', '\Sigma^{5} Y')
 
%end

%for kplot = 1:4
%  subplot(2,2,kplot)

 %       [N001, histbin001] = hist(pwr001(:,binno,kplot));
%        [N002, histbin002] = hist(pwr002(:,binno,kplot));
%        [N005, histbin005] = hist(pwr005(:,binno,kplot));
%        [N002sum, histbin002sum] = hist(pwr002sum(:,binno,kplot));
%        [N005sum, histbin005sum] = hist(pwr005sum(:,binno,kplot));
        
%        totalarea001 = trapz(histbin001, N001);
%        totalarea002 = trapz(histbin002, N002);
%        totalarea005 = trapz(histbin005, N005);
%        totalarea_sum002 = trapz(histbin002sum, N002sum);
%        totalarea_sum005 = trapz(histbin005sum, N005sum);
        
%        plot(histbin001, log(N001/(totalarea001)), '--k', ...
%             histbin002, log(N002/(totalarea002)), '-k', ...
%             histbin005, log(N005/(totalarea005)), ':k', ...
%             histbin002sum, log(N002sum/(totalarea_sum002)), 'ok', ...
%             histbin005sum, log(N005sum/(totalarea_sum005)), 'xk')
%ylabel('ln p(Y)')
%  legend('T_I = 1 ms', 'T_I = 2 ms', 'T_I = 5 ms', '\Sigma^{2} Y', '\Sigma^{5} Y')
  
%end

%
% Generate cross-correlations and correlatio times - copied from corrtime.m
% script. 
%

ndelays = size(delsyn_integ001t,2);
ncases = size(wfsyn_integ001t,3);

C001 = zeros(21,ndelays);
C002 = zeros(21,ndelays);
C003 = zeros(21,ndelays);
C005 = zeros(21,ndelays);
C01 = zeros(21,ndelays);
C001s = zeros(21,ndelays);
C002s = zeros(21,ndelays);
C003s = zeros(21,ndelays);
C005s = zeros(21,ndelays);
C01s = zeros(21,ndelays);

lags001 = zeros(201,ndelays);
lags002 = zeros(201,ndelays);
lags005 = zeros(201,ndelays);
lags01 = zeros(201,ndelays);
lags001s = zeros(201,ndelays);
lags002s = zeros(201,ndelays);
lags005s = zeros(201,ndelays);
lags01s = zeros(201,ndelays);

corrtime001 = zeros(ndelays,ncases);
corrtime002 = zeros(ndelays,ncases);
corrtime003 = zeros(ndelays,ncases);
corrtime005 = zeros(ndelays,ncases);
corrtime01 = zeros(ndelays,ncases);
corrtime001s = zeros(ndelays,ncases);
corrtime002s = zeros(ndelays,ncases);
corrtime003s = zeros(ndelays,ncases);
corrtime005s = zeros(ndelays,ncases);
corrtime01s = zeros(ndelays,ncases);

corrtime001m = zeros(ndelays,ncases);
corrtime002m = zeros(ndelays,ncases);
corrtime003m = zeros(ndelays,ncases);
corrtime005m = zeros(ndelays,ncases);
corrtime01m = zeros(ndelays,ncases);
corrtime001sm = zeros(ndelays,ncases);
corrtime002sm = zeros(ndelays,ncases);
corrtime003sm = zeros(ndelays,ncases);
corrtime005sm = zeros(ndelays,ncases);
corrtime01sm = zeros(ndelays,ncases);

corrtime001_fit = zeros(ndelays,ncases);
corrtime002_fit = zeros(ndelays,ncases);
corrtime003_fit = zeros(ndelays,ncases);
corrtime005_fit = zeros(ndelays,ncases);
corrtime01_fit = zeros(ndelays,ncases);
corrtime001s_fit = zeros(ndelays,ncases);
corrtime002s_fit = zeros(ndelays,ncases);
corrtime003s_fit = zeros(ndelays,ncases);
corrtime005s_fit = zeros(ndelays,ncases);
corrtime01s_fit = zeros(ndelays,ncases);

corrtime001m_fit = zeros(ndelays,ncases);
corrtime002m_fit = zeros(ndelays,ncases);
corrtime003m_fit = zeros(ndelays,ncases);
corrtime005m_fit = zeros(ndelays,ncases);
corrtime01m_fit = zeros(ndelays,ncases);
corrtime001sm_fit = zeros(ndelays,ncases);
corrtime002sm_fit = zeros(ndelays,ncases);
corrtime003sm_fit = zeros(ndelays,ncases);
corrtime005sm_fit = zeros(ndelays,ncases);
corrtime01sm_fit = zeros(ndelays,ncases);


resnorm001m = zeros(ndelays,ncases);
resnorm002m = zeros(ndelays,ncases);
resnorm003m = zeros(ndelays,ncases);
resnorm005m = zeros(ndelays,ncases);

resnorm001sm = zeros(ndelays,ncases);
resnorm002sm = zeros(ndelays,ncases);
resnorm003sm = zeros(ndelays,ncases);
resnorm005sm = zeros(ndelays,ncases);

residual001m = zeros(ndelays,ncases);
residual002m = zeros(ndelays,ncases);
residual003m = zeros(ndelays,ncases);
residual005m = zeros(ndelays,ncases);

residual001sm = zeros(ndelays,ncases);
residual002sm = zeros(ndelays,ncases);
residual003sm = zeros(ndelays,ncases);
residual005sm = zeros(ndelays,ncases);


for kcase = 1:ncases
   for k=1:ndelays
       
   fprintf(' Case = %4i   Delay = %4i \n', kcase, k)
  % [C001s(:,k), lags001s] = xcorr(squeeze(wfsum001_short(:,k,2)), 100,'unbiased');
   [C001(:,k), lags001] = xcorr(squeeze(wfsyn_integ001t(:,k,kcase)), 10, 'unbiased');
   [C002(:,k), lags002] = xcorr(squeeze(wfsyn_integ002t(:,k,kcase)), 10, 'unbiased');
   [C002s(:,k), lags002s] = xcorr(squeeze(wfsum002(:,k,kcase)), 10, 'unbiased');
   [C003(:,k), lags003] = xcorr(squeeze(wfsyn_integ003t(:,k,kcase)), 10, 'unbiased');
   [C003s(:,k), lags003s] = xcorr(squeeze(wfsum003(:,k,kcase)), 10, 'unbiased');
   [C005s(:,k), lags005s] = xcorr(squeeze(wfsum005(:,k,kcase)), 10, 'unbiased');
   [C005(:,k), lags005] = xcorr(squeeze(wfsyn_integ005t(:,k,kcase)), 10, 'unbiased');
   [C01s(:,k), lags01s] = xcorr(squeeze(wfsum01(:,k,kcase)), 10, 'unbiased');
   [C01(:,k), lags01] = xcorr(squeeze(wfsyn_integ01t(:,k,kcase)), 10, 'unbiased');
   end
   corrtime001(:,kcase) = trapz(lags001', real(C001))./max(real(C001),[],1);
   corrtime002(:,kcase) = trapz(lags002'*2, real(C002))./max(real(C002),[],1);
   corrtime002s(:,kcase) = trapz(lags002s'*2, real(C002s))./max(real(C002s),[],1);
   corrtime003(:,kcase) = trapz(lags003'*3, real(C003))./max(real(C003),[],1);
   corrtime003s(:,kcase) = trapz(lags003s'*3, real(C003s))./max(real(C003s),[],1);
   corrtime005(:,kcase) = trapz(lags005'*5, real(C005))./max(real(C005),[],1);
   corrtime005s(:,kcase) = trapz(lags005s'*5, real(C005s))./max(real(C005s),[],1);
   
   corrtime001_fit(:,kcase) = chartime(lags001, real(C001));
   corrtime002_fit(:,kcase) = chartime(lags002*2, real(C002));
   corrtime002s_fit(:,kcase) = chartime(lags002s*2, real(C002s));
   corrtime003_fit(:,kcase) = chartime(lags003*3, real(C003));
   corrtime003s_fit(:,kcase) = chartime(lags003s*3, real(C003s));
   corrtime005_fit(:,kcase) = chartime(lags005*5, real(C005));
   corrtime005s_fit(:,kcase) = chartime(lags005s*5, real(C005s));
   
   corrtime001m(:,kcase) = trapz(lags001', abs(C001))./max(abs(C001),[],1);
   corrtime002m(:,kcase) = trapz(lags002'*2, abs(C002))./max(abs(C002),[],1);
   corrtime002sm(:,kcase) = trapz(lags002s'*2, abs(C002s))./max(abs(C002s),[],1);
   corrtime003m(:,kcase) = trapz(lags003'*3, abs(C003))./max(abs(C003),[],1);
   corrtime003sm(:,kcase) = trapz(lags003s'*3, abs(C003s))./max(abs(C003s),[],1);
   corrtime005m(:,kcase) = trapz(lags005'*5, abs(C005))./max(abs(C005),[],1);
   corrtime005sm(:,kcase) = trapz(lags005s'*5, abs(C005s))./max(abs(C005s),[],1);
   
   [corrtime001m_fit(:,kcase), resnorm001m(:,kcase), res001m(:,kcase)] = chartime(lags001, abs(C001));
   [corrtime002m_fit(:,kcase),resnorm002m(:,kcase), res002m(:,kcase)] = chartime(lags002*2, abs(C002));
   [corrtime002sm_fit(:,kcase), resnorm002sm(:,kcase), res002sm(:,kcase)] = chartime(lags002s*2, abs(C002s));
   [corrtime003m_fit(:,kcase), resnorm003m(:,kcase), res003m(:,kcase)] = chartime(lags003*3, abs(C003));
   [corrtime003sm_fit(:,kcase), resnorm003sm(:,kcase), res003sm(:,kcase)] = chartime(lags003s*3, abs(C003s));
   [corrtime005m_fit(:,kcase), resnorm005m(:,kcase), res005m(:,kcase)] = chartime(lags005*5, abs(C005));
   [corrtime005sm_fit(:,kcase), resnorm005sm(:,kcase), res005sm(:,kcase)] = chartime(lags005s*5, abs(C005s));
   
end

subplotid = {'(a)', '(b)', '(c)', '(d)'};
%set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

% Change default text fonts.
%set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 12)

axisset = [-1.5 2.5 0 4;
           -1.5 2.5 0 4;
           -1.5 2.5 0 5; 
           -1.5 2.5 0 7];

figure(3)

for kplot=1:4
  subplot(2,2,kplot)
  plot(delsyn_integ001t(1,:,kplot), corrtime001m_fit(:,kplot), ':k', ...
  delsyn_integ002t(1,:,kplot), corrtime002m_fit(:,kplot), '--k', ...
  delsyn_integ002t(1,:,kplot), corrtime002sm_fit(:,kplot), 'ok', ...
  delsyn_integ003t(1,:,kplot), corrtime003m_fit(:,kplot), '-k', ...
  delsyn_integ003t(1,:,kplot), corrtime003sm_fit(:,kplot), 'xk')
  xlabel('Delay (chips)')
  ylabel('Correlation Time (ms)')
  axis(axisset(kplot,:))
  title(subplotid(kplot))
end

legend('T_I = 1 ms', 'T_I = 2 ms', 'T_I = 1 ms, N_C=2', 'T_I = 3 ms', 'T_I = 1 ms, N_C=3')

%figure(4)

%for kplot=1:4
%  subplot(2,2,kplot)
%  plot(delsyn_integ001t(1,:,kplot), corrtime001_fit(:,kplot), ':k', ...
%  delsyn_integ002t(1,:,kplot), corrtime002_fit(:,kplot), '--k', ...
%  delsyn_integ002t(1,:,kplot), corrtime002s_fit(:,kplot), 'o', ...
%  delsyn_integ003t(1,:,kplot), corrtime003_fit(:,kplot), '-k', ...
%  delsyn_integ003t(1,:,kplot), corrtime003s_fit(:,kplot), 'x')
%  xlabel('Delay (chips)')
%  ylabel('Correlation Time (ms)')
%    axis([-1.7 2.7 0 5])
%  title(subplotid(kplot))
%end

%legend('T_I = 1 ms', 'T_I = 2 ms', 'T_I = 1 ms, N_C=2', 'T_I = 3 ms', 'T_I = 1 ms, N_C=3')

figure(5)

for kplot=1:4
  subplot(2,2,kplot)
  plot(delsyn_integ001t(1,:,kplot), resnorm001m(:,kplot), ':k', ...
  delsyn_integ002t(1,:,kplot), resnorm002m(:,kplot), '--k', ...
  delsyn_integ002t(1,:,kplot), resnorm002sm(:,kplot), 'ok', ...
  delsyn_integ003t(1,:,kplot), resnorm003m(:,kplot), '-k', ...
  delsyn_integ003t(1,:,kplot), resnorm003sm(:,kplot), 'xk')
  xlabel('Delay (chips)')
  ylabel('Residual norm')
  %  axis([-1.7 2.7 0 5])
  title(subplotid(kplot))
  
end

legend('T_I = 1 ms', 'T_I = 2 ms', 'T_I = 1 ms, N=2', 'T_I = 3 ms', 'T_I = 1 ms, N=3')