%
% Script to generate a plot of the spectra, and the derived roughness
% models for comparison, and test vs.  what is in the papers
%
%  Dec 18, 2007 at Starlab
%
%  Elfouhaily spectrum
%
wavenum = power(10,-3:0.01:4);
omega = 0.84;
%U10 = 3:2:21;  %Same numbers as in the paper
U10 = 3:3:21;  % smaller range for plotting in book. 
nws = size(U10,2);

Smat = zeros(nws,size(wavenum,2));

for k=1:nws
  Smat(k,:) = elfouhaily1d(wavenum, U10(k), omega);
end


subplot(1,1,1)
loglog(wavenum, Smat,[11   11], [1e-16, 1e3] )
axis([1e-3 1e4 1e-15 1e3])
ha1 = gca;
xlabel('Wavenumber k (rad/m)')
ylabel('Omnidirectional Elevation Spectrum: S(k)')

%title('These plots should match Figure 8')
grid

pause

subplot(1,1,1)
loglog(wavenum, Smat .* (ones(nws,1)*wavenum).^3 );
xlabel('Wavenumber k (rd/m)')
ylabel('Curvature Spectrum B(k)=k^3 S(k)')
axis([1e-3 1e4 1e-4 1])
ha2 = gca;
title('These plots should match Figure 8')
grid

pause

for k=1:nws
  [sig2u(k,:), sig2c(k,:), sig2o(k,:)] = var_from_spec(U10(k), omega, 1e3);
end

gcpcm = coxmunk_pdf( U10', 40);
gcpwu = coxmunk_pdf( U10', 40, 'wu');

subplot(1,1,1)
plot(U10, gcpcm(:,1)+gcpcm(:,2), '*', U10, gcpwu(:,1)+gcpwu(:,2), 'x', ...
     U10, sig2u + sig2c, U10, sig2o, '--')
xlabel('U_{10} (m/s)')
ylabel('MSS')
axis([0 24 0 0.12])
legend('Cox & Munk (1957)', 'Wu (1992)', 'Integration of k^2S(k)')
title('This plot sould match figure 7')