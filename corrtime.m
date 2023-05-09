%
% Generate cross-correlations and correlatio times.
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
   
   corrtime001m_fit(:,kcase) = chartime(lags001, abs(C001));
   corrtime002m_fit(:,kcase) = chartime(lags002*2, abs(C002));
   corrtime002sm_fit(:,kcase) = chartime(lags002s*2, abs(C002s));
   corrtime003m_fit(:,kcase) = chartime(lags003*3, abs(C003));
   corrtime003sm_fit(:,kcase) = chartime(lags003s*3, abs(C003s));
   corrtime005m_fit(:,kcase) = chartime(lags005*5, abs(C005));
   corrtime005sm_fit(:,kcase) = chartime(lags005s*5, abs(C005s));
   
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

figure(1)

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

legend('T_I = 1 ms', 'T_I = 2 ms', 'T_I = 1 ms, N=2', 'T_I = 3 ms', 'T_I = 1 ms, N=3')

figure(2)

for kplot=1:4
  subplot(2,2,kplot)
  plot(delsyn_integ001t(1,:,kplot), corrtime001_fit(:,kplot), ':k', ...
  delsyn_integ002t(1,:,kplot), corrtime002_fit(:,kplot), '--k', ...
  delsyn_integ002t(1,:,kplot), corrtime002s_fit(:,kplot), 'o', ...
  delsyn_integ003t(1,:,kplot), corrtime003_fit(:,kplot), '-k', ...
  delsyn_integ003t(1,:,kplot), corrtime003s_fit(:,kplot), 'x')
  xlabel('Delay (chips)')
  ylabel('Correlation Time (ms)')
    axis([-1.7 2.7 0 5])
  title(subplotid(kplot))
end

legend('T_I = 1 ms', 'T_I = 2 ms', 'T_I = 1 ms, N=2', 'T_I = 3 ms', 'T_I = 1 ms, N=3')