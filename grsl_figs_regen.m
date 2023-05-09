%
% Regenerate the figures for the GRSL paper.
%
% Modified 2/14 by JLG to generate binned statistics for simulator paper
%
%
% Everything until "new part" is unchanged and should produce the same
% figures as in the GRSL paper
%
clear all

load ./2Mar09_L2DATA_final_jlg

offaxismask = logical(cosgam < 0.851);

%
% Figure 3 - example waveform
% 
figure(3)

plot(del_inp(450*300:451*300,:), wf_inp(450*300:451*300,:), '.k', ...
    chips_model(451,:), wf_model(451,:), 'k')
xlabel('Delay: \tau (chips)')
ylabel('Post-correlation waveform: |Y(\tau, f_{d0})|^2')
grid
axis([-2 6 -0.1 0.7])

2*xest(451,1)^2
cosgam(451)

%
% Run the empirical model derived from J Voo
%

%tbvfit = fitevalfvb(2*xest(~offaxismask,1).^2);
% tbhfit = fitevalfhb(2*xest(~offaxismask,1).^2);

%
% Run the emprical model derived by JLG 12/2010
%
x2plot = 0.015:0.001:0.037;

[tbvplot, tbhplot] = poly4model( x2plot);
[tbvfit, tbhfit] = poly4model( 2*xest(~offaxismask,1).^2);


%
% Figure 4 - TB vs. MSS
%
figure(4)

subplot(2,1,1)
plot(2* xest(~offaxismask,1).^2, tbv(~offaxismask), 'ok', ...
    2* xest(offaxismask,1).^2, tbv(offaxismask), 'xk', ...
        x2plot, tbvplot, 'k')
axis([0.015 0.04 120 128])
grid
xlabel('\sigma_{iso}^2')
ylabel('T_{bv} (K)')


subplot(2,1,2)
plot(2* xest(~offaxismask,1).^2, tbh(~offaxismask), 'ok', ...=-
    2* xest(offaxismask,1).^2, tbh(offaxismask), 'xk', ...
           x2plot, tbhplot, 'k')
axis([0.015 0.04 74 84])
grid
xlabel('\sigma_{iso}^2')
ylabel('T_{bh} (K)')

%
% Figure 5 - scatterplot
%
figure(5)

subplot(1,1,1)
plot(tbv(~offaxismask), tbvfit,'.k', [120 128], [120 128], 'k')
axis([120 128 120 128])
xlabel('PALS    T_{bv} (K)')
ylabel('GNSS-R    T_{bv} (K)')
grid
axis square

figure(6) 

subplot(1,1,1)
plot(tbh(~offaxismask), tbhfit,  '.k', [74 83], [74 83], 'k')
axis([73 83 73 83])
xlabel('PALS    T_{bh} (K)')
ylabel('GNSS-R    T_{bh} (K)')
grid
axis square

stdresv = std(tbvfit - tbv(~offaxismask))

stdresh = std(tbhfit - tbh(~offaxismask))

meanresv = mean(tbvfit - tbv(~offaxismask))

meanresh = mean(tbhfit - tbh(~offaxismask))

rmsresv = sqrt(mean((tbvfit - tbv(~offaxismask)).^2))

rmsresh = sqrt(mean((tbhfit - tbh(~offaxismask)).^2))

%
%  New part - Bin data. 
% 
tbv_pals = tbv(~offaxismask);
tbh_pals = tbh(~offaxismask);

tbv_error = tbvfit - tbv_pals;
tbh_error = tbhfit - tbh_pals; 

figure(7)
subplot(2,1,1)
plot(tbv_pals, tbv_error, 'x')
subplot(2,1,2)
plot(tbh_pals, tbh_error, 'x')



