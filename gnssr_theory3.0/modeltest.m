%
% Model consistency test for GNSS-R waveform generator 
%
clear
load mpcp
fprintf(' Model consistency test ... \n')
fprintf(' Receiver altitude = %12.5f m \n', cp.alt_m);
fprintf(' SV elevation = %12.5f deg \n', cp.gammadeg);
fprintf(' Receiver velocity = %12.5f, %12.5f, %12.5f m/s \n', cp.VR(1), cp.VR(2), cp.VR(3));
fprintf(' GPS velocity = %12.5f, %12.5f, %12.5f m/s \n', cp.VG(1), cp.VG(2), cp.VG(3));

dftilde = mp.dftilde; 
ftilde_max = mp.ftilde_max; 
dttilde = mp.dttilde;
ttilde_max = mp.ttilde_max; 

samplelag1 = 100;  % define 4 points in the waveform where spectrum and acf are compared
samplelag2 = 180;
samplelag3 = 250;
samplelag4 = 400;


nwf_test = 20000;  % number of synthetic waveform to generate to test. 

fprintf(' Delta f-tilde =  %12.5f Hz \n', dftilde);
fprintf(' f-tilde range = +/- %12.5f Hz \n', ftilde_max);

fprintf(' Delta t-tilde =  %12.5f ms \n', dttilde*1000);
fprintf(' t-tilde range = +/- %12.5f ms \n', ttilde_max*1000);

fprintf(' Delta theta = %12.5 Hz \n', mp.thetastep);

% 
% Compute the spectra for the time series from each bin.
%
fprintf(' Computing the power spectra ... \n');


ftilde = [-ftilde_max:dftilde:ftilde_max ];  % Hz
nftilde = size(ftilde,2);

%fdsurf = zeros(1334,629,nftilde);

tic

for k=1:size(ftilde,2)
 [dummy, pcdarray, corrspec(:,k), pwr, fs, fdsurf ] = ...
     wf_from_gcp_dopp( 0, cp, mp, 0, ftilde(k), dftilde);
 fprintf('   ftilde = %12.5f Hz \n', ftilde(k))
end

fd_runtime = toc;

tic
[wf_array, pcdarray, dummy, pwr, fs, fdsurf ] = ...
     wf_from_gcp_dopp( 0, cp, mp, 0, 0, 0);
 wf_runtime = toc;

spec_int = trapz(ftilde, corrspec,2);


%
% Compute the autocorrelation for the time series from each bin.
%
fprintf(' Computing the autocorrelation  ... \n');


ttilde = [-ttilde_max:dttilde:ttilde_max];  % sec
nttilde = size(ttilde,2);

%fdsurf = zeros(1334,629,nftilde);

tic

for k=1:size(ttilde,2)
 [yy(:,k), pcdarray, dummy, pwr, fs, fdsurf ] = ...
     wf_from_gcp_dopp( 0, cp, mp, ttilde(k), 0, dftilde);
 fprintf('   ttilde = %12.5f ms \n', ttilde(k)*1000)
end

td_runtime  = toc;

subplot(1,1,1)
plot(pcdarray, wf_array, pcdarray, spec_int, pcdarray, yy(:,(nttilde+1)/2));
xlabel('Delay \tau (m)')
ylabel('Waveform (unity input power)')
legend('Y^2(\tau)','\int W_Y(f) df', 'R_Y(0)')
pause

fprintf(' Test that R_Y(0) is real:  %4i \n', isreal(yy(:,(nttilde+1)/2)));

fprintf(' Test consistency between time- and frequency-domain models \n')

%
% FFT the autocorrelation
%
ftilde_fft = -(1/(2*dttilde)):(1/(2*ttilde_max)):(1/(2*dttilde));
corrspec_fft = fftshift(fft(ifftshift(yy),[],2))*dttilde;

subplot(1,1,1)
plot(ftilde, corrspec(200,:), ftilde_fft, real(corrspec_fft(200,:)));
xlabel('f (Hz)')
ylabel('W_Y(f)')
legend('W_Y(f) model', 'FFT(R_Y(\tau)) model')

fprintf(' Runtimes: \n')
fprintf('  W_Y(f): %12.5f sec\n', fd_runtime);
fprintf('  R_Y(tau): %12.5f sec\n', td_runtime);
fprintf('  Y(tau): %12.5f sec\n', wf_runtime);

fprintf('\n Generating synthetic data ... \n')
fprintf('   %8i wavforms simulated \n', nwf_test);
fprintf('   No bin-bin correlation \n')
fprintf('   C/N0 = %f12.5 dB \n', cp.CN0dB); 

[cdsim, y, dummy1, dummy2, dummy3] = wfsim(mp, cp, power(10,cp.CN0dB/10), nwf_test, 0, 0);

nf = mean(y(:,1).*conj(y(:,1)));

subplot(1,1,1)
plot(pcdarray, wf_array, cdsim, mean(y.*conj(y)-nf,1), ...
    [pcdarray(samplelag1), pcdarray(samplelag2), pcdarray(samplelag3), pcdarray(samplelag4)], ...
    [wf_array(samplelag1), wf_array(samplelag2), wf_array(samplelag3), wf_array(samplelag4)], 'o')
xlabel('Delay \tau (m)')
ylabel(' Y^2(\tau)')
legend('Numerical Model','Simulated Data','Sample points for Spectrum comparison')

pause

fprintf(' Computing spectrum ... \n')

BW = 1/cp.Ti;  % sample rate, also bandwidth of noise.

[wysim1, ftildesim1] = pwelch(y(:,samplelag1), [], [], [], BW);
[wysim2, ftildesim2] = pwelch(y(:,samplelag2), [], [], [], BW);
[wysim3, ftildesim3] = pwelch(y(:,samplelag3), [], [], [], BW);
[wysim4, ftildesim4] = pwelch(y(:,samplelag4), [], [], [], BW);

nspec = size(ftildesim1,1);
nfspec1 = mean(wysim1(nspec/2-256:nspec/2+256));  % compute noise at mid-point (not "flipped")
nspec = size(ftildesim2,1);
nfspec2 = mean(wysim2(nspec/2-256:nspec/2+256));  % compute noise at mid-point (not "flipped")
nspec = size(ftildesim3,1);
nfspec3 = mean(wysim3(nspec/2-256:nspec/2+256));  % compute noise at mid-point (not "flipped")
nspec = size(ftildesim4,1);
nfspec4 = mean(wysim4(nspec/2-256:nspec/2+256));  % compute noise at mid-point (not "flipped")

fprintf('Noise floor = %12.5f \n', nf)
fprintf('White noise spectrum = (nf / BW) = %12.8f \n', nf/BW);
fprintf(' For lag number 1: %12.8f \n', nfspec1);
fprintf(' For lag number 2: %12.8f \n', nfspec2);
fprintf(' For lag number 3: %12.8f \n', nfspec3);
fprintf(' For lag number 4: %12.8f \n', nfspec4);


subplot(2,2,1)
plot(ftildesim1-max(ftildesim1)/2, fftshift(wysim1)-nfspec1, '.y', ftilde, corrspec(samplelag1,:),'-x')
title(['Lag number:', num2str(samplelag1), ' Delay = ', num2str(pcdarray(samplelag1)),' m '])
xlabel('f-tilde (Hz)')
ylabel('W_Y(f)')
legend('Simulated','Model')

subplot(2,2,2)
plot(ftildesim2-max(ftildesim2)/2, fftshift(wysim2)-nfspec2, '.y', ftilde, corrspec(samplelag2,:),'-x')
title(['Lag number:', num2str(samplelag2), ' Delay = ', num2str(pcdarray(samplelag2)),' m '])
xlabel('f-tilde (Hz)')
ylabel('W_Y(f)')

subplot(2,2,3)
plot(ftildesim3-max(ftildesim3)/2, fftshift(wysim3)-nfspec3, '.y', ftilde, corrspec(samplelag3,:),'-x')
title(['Lag number:', num2str(samplelag3), ' Delay = ', num2str(pcdarray(samplelag3)),' m '])
xlabel('f-tilde (Hz)')
ylabel('W_Y(f)')

subplot(2,2,4)
plot(ftildesim4-max(ftildesim4)/2, fftshift(wysim4)-nfspec4, '.y', ftilde, corrspec(samplelag4,:),'-x')
title(['Lag number:', num2str(samplelag4), ' Delay = ', num2str(pcdarray(samplelag4)),' m '])
xlabel('f-tilde (Hz)')
ylabel('W_Y(f)')



fprintf(' Computing the autocorrelation ... \n')

[Rysim1, lags1] = xcorr(y(:,samplelag1), y(:,samplelag1), 100, 'unbiased');
[Rysim2, lags2] = xcorr(y(:,samplelag2), y(:,samplelag2), 100, 'unbiased');
[Rysim3, lags3] = xcorr(y(:,samplelag3), y(:,samplelag3), 100, 'unbiased');
[Rysim4, lags4] = xcorr(y(:,samplelag4), y(:,samplelag4), 100, 'unbiased');

Rymidpoint = (size(Rysim1,1) +1)/2;

% Remove noise floor at 0-lag:

Rysim1(Rymidpoint) = Rysim1(Rymidpoint) -nf;
Rysim2(Rymidpoint) = Rysim2(Rymidpoint) -nf;
Rysim3(Rymidpoint) = Rysim3(Rymidpoint) -nf;
Rysim4(Rymidpoint) = Rysim4(Rymidpoint) -nf;


subplot(2,2,1)
plot(lags1*cp.Ti, real(Rysim1), '-x', ttilde, real(yy(samplelag1,:)), '-') 
title(['Lag number:', num2str(samplelag1), ' Delay = ', num2str(pcdarray(samplelag1)),' m '])
xlabel('t-tilde (ms)')
ylabel('Real[R_Y(t)]')
legend('Simulated','Model')
grid

subplot(2,2,2)
plot(lags2*cp.Ti, real(Rysim2), '-x', ttilde, real(yy(samplelag2,:)), '-') 
title(['Lag number:', num2str(samplelag2), ' Delay = ', num2str(pcdarray(samplelag2)),' m '])
xlabel('t-tilde (sec)')
ylabel('Real[R_Y(t)]')
grid


subplot(2,2,3)
plot(lags3*cp.Ti, real(Rysim3), '-x', ttilde, real(yy(samplelag3,:)), '-') 
title(['Lag number:', num2str(samplelag3), ' Delay = ', num2str(pcdarray(samplelag3)),' m '])
xlabel('t-tilde (sec)')
ylabel('Real[R_Y(t)]')
grid

subplot(2,2,4)
plot(lags4*cp.Ti, real(Rysim4), '-x', ttilde, real(yy(samplelag4,:)), '-') 
title(['Lag number:', num2str(samplelag4), ' Delay = ', num2str(pcdarray(samplelag4)),' m '])
xlabel('t-tilde (sec)')
ylabel('Real[R_Y(t)]')
grid

pause

subplot(2,2,1)
plot(lags1*cp.Ti, imag(Rysim1), '-x', ttilde, imag(yy(samplelag1,:)), '-') 
title(['Lag number:', num2str(samplelag1), ' Delay = ', num2str(pcdarray(samplelag1)),' m '])
xlabel('t-tilde (ms)')
ylabel('imag[R_Y(t)]')
legend('Simulated','Model')
grid

subplot(2,2,2)
plot(lags2*cp.Ti, imag(Rysim2), '-x', ttilde, imag(yy(samplelag2,:)), '-') 
title(['Lag number:', num2str(samplelag2), ' Delay = ', num2str(pcdarray(samplelag2)),' m '])
xlabel('t-tilde (sec)')
ylabel('imag[R_Y(t)]')
grid

subplot(2,2,3)
plot(lags3*cp.Ti, imag(Rysim3), '-x', ttilde, imag(yy(samplelag3,:)), '-') 
title(['Lag number:', num2str(samplelag3), ' Delay = ', num2str(pcdarray(samplelag3)),' m '])
xlabel('t-tilde (sec)')
ylabel('imag[R_Y(t)]')
grid

subplot(2,2,4)
plot(lags4*cp.Ti, imag(Rysim4), '-x', ttilde, imag(yy(samplelag4,:)), '-') 
title(['Lag number:', num2str(samplelag4), ' Delay = ', num2str(pcdarray(samplelag4)),' m '])
xlabel('t-tilde (sec)')
ylabel('imag[R_Y(t)]')
grid

pause

subplot(2,2,1)
plot(lags1*cp.Ti, abs(Rysim1), '-x', ttilde, abs(yy(samplelag1,:)), '-') 
title(['Lag number:', num2str(samplelag1), ' Delay = ', num2str(pcdarray(samplelag1)),' m '])
xlabel('t-tilde (ms)')
ylabel('abs[R_Y(t)]')
legend('Simulated','Model')
grid

subplot(2,2,2)
plot(lags2*cp.Ti, abs(Rysim2), '-x', ttilde, abs(yy(samplelag2,:)), '-') 
title(['Lag number:', num2str(samplelag2), ' Delay = ', num2str(pcdarray(samplelag2)),' m '])
xlabel('t-tilde (sec)')
ylabel('abs[R_Y(t)]')
grid

subplot(2,2,3)
plot(lags3*cp.Ti, abs(Rysim3), '-x', ttilde, abs(yy(samplelag3,:)), '-') 
title(['Lag number:', num2str(samplelag3), ' Delay = ', num2str(pcdarray(samplelag3)),' m '])
xlabel('t-tilde (sec)')
ylabel('abs[R_Y(t)]')
grid

subplot(2,2,4)
plot(lags4*cp.Ti, abs(Rysim4), '-x', ttilde, abs(yy(samplelag4,:)), '-') 
title(['Lag number:', num2str(samplelag4), ' Delay = ', num2str(pcdarray(samplelag4)),' m '])
xlabel('t-tilde (sec)')
ylabel('abs[R_Y(t)]')
grid

pause

subplot(2,2,1)
plot(lags1*cp.Ti, angle(Rysim1), '-x', ttilde, angle(yy(samplelag1,:)), '-') 
title(['Lag number:', num2str(samplelag1), ' Delay = ', num2str(pcdarray(samplelag1)),' m '])
xlabel('t-tilde (ms)')
ylabel('angle[R_Y(t)]')
legend('Simulated','Model')
grid

subplot(2,2,2)
plot(lags2*cp.Ti, angle(Rysim2), '-x', ttilde, angle(yy(samplelag2,:)), '-') 
title(['Lag number:', num2str(samplelag2), ' Delay = ', num2str(pcdarray(samplelag2)),' m '])
xlabel('t-tilde (sec)')
ylabel('angle[R_Y(t)]')
grid


subplot(2,2,3)
plot(lags3*cp.Ti, angle(Rysim3), '-x', ttilde, angle(yy(samplelag3,:)), '-') 
title(['Lag number:', num2str(samplelag3), ' Delay = ', num2str(pcdarray(samplelag3)),' m '])
xlabel('t-tilde (sec)')
ylabel('angle[R_Y(t)]')
grid

subplot(2,2,4)
plot(lags4*cp.Ti, angle(Rysim4), '-x', ttilde, angle(yy(samplelag4,:)), '-') 
title(['Lag number:', num2str(samplelag4), ' Delay = ', num2str(pcdarray(samplelag4)),' m '])
xlabel('t-tilde (sec)')
ylabel('angle[R_Y(t)]')
grid


%[wf_array, pcdarray, dummy, pwr, fs, fdsurf, Rtau] = ...
% wf_from_gcp_dopp( 0, cp, mp, 0, 0, dftilde);


return
