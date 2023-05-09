%
% Generate an illustration of the change of variables to delay-Doppler
%
% This will run MATLAB to generate a basic illustration which relaistic
% number, this can be modified in a drawing progpwdram to clean it up
%
% Will generate the example for case "caseno" of the simulator run
%
clear

load 2Mar09_level2
latexfiledir = '/Users/jgarriso/jgarriso/pubs/papers/2012_tgars_simulator/draft3/'

caseno = 1;
dthetadeg = 5;  % deg
dtheta = dthetadeg*pi/180;  % radian
ddelay = 30; % meters
ftilde = 100; % Hz.
delaydopp_tol = 100; % meters - how close delay-Doppler intersction needs to be to plot a point

elevation = elhist(caseno);  % deg


altitude = althist(caseno);% meters
VR =  VRscat(caseno,:);
VG = VGscat(caseno,:);

VR(3) = 0;

delay1chips = 0.75; % chips
delay1 = delay1chips*293; % meters
delay2 = delay1 + ddelay; % small delta delay equial to one step in the integrator
delaychips = [0.25 0.5 0.75 1.0 1.25 1.5]; % chips
delays = delaychips* 293;  % meters

dopp1 = 0; % Hz
dopp2 = 200;  %Hz
dopps = -500:50:500;  %Hz
ndopp = size(dopps,2);

thetas = [-180:30:180]*pi/180;
theta1deg = 30;  % deg
theta1 = theta1deg*pi/180;  % rad
theta2 = theta1 + dtheta;

rmaxtheta = 6000;
nthetas = size(thetas,1);

[x1, y1, r, theta_isorange1, x10, xs, a, b] = isorange( altitude, elevation, delay1);
[x2, y2, r, theta_dummy, x20, xs, a, b] = isorange( altitude, elevation, delay2);
[x, y, r, theta, x20, xs, a, b] = isorange( altitude, elevation, delays);

[xd1, yd1, thetad1, betad] = ...
         isodopp( altitude, elevation, dopp1, VG, VR, pi-0.01); 
[xd2, yd2, thetad2, betad] = ...
         isodopp( altitude, elevation, dopp2, VG, VR, pi-0.01); 
%[xd1, yd1, thetad, betad] = ...
% isodopp( altitude, elevation, dopp1, VG, VR, pi-0.01);     
[xd, yd, thetad, betad,fs, gotsoln] = ...
         isodopp( altitude, elevation, dopps, VG, VR, pi-0.01);

%
% Location of Doppler intersections on iso-range ellipse
%
[doppoff, doppindex] = min( abs(dopps-ftilde));
xd_example = xd(gotsoln(:,doppindex),doppindex);
yd_example = yd(gotsoln(:,doppindex),doppindex);

delaydopp_vecx = ones(size(xd_example,1),1)*x1' - xd_example * ones(1,size(x1,1));
delaydopp_vecy = ones(size(yd_example,1),1)*y1' - yd_example * ones(1,size(y1,1));

[delaydopp_nearest, delaydopp_nearest_index] = min( delaydopp_vecx.^2+delaydopp_vecy.^2, [], 1);

delaydopp_intersect = delaydopp_nearest < delaydopp_tol;
delaydopp_theta_intersect = thetad(delaydopp_nearest_index(delaydopp_intersect));

delaydopp_intersect_x = x1(delaydopp_intersect);
delaydopp_intersect_y =y1(delaydopp_intersect);

     
xtheta = rmaxtheta * cos(thetas);
ytheta = rmaxtheta * sin(thetas);

xtheta1 = rmaxtheta * cos(theta1);
ytheta1 = rmaxtheta * sin(theta1);
xtheta2 = rmaxtheta * cos(theta2);
ytheta2 = rmaxtheta * sin(theta2);

xpolar = [x10*ones(1,13); xtheta+x10];
ypolar = [0*ones(1,13); ytheta];

xpolar1 = [x10; xtheta1+x10];
ypolar1 = [0; ytheta1];
xpolar2 = [x10; xtheta2+x10];
ypolar2 = [0; ytheta2];

%plot(x1, y1, 'k', x10, 0, 'ok', xs, 0, '*k', 0, 0, 'xk', ...
%     x2, y2, 'r', x20, 0, 'or', xs, 0, '*k', 0, 0, 'xk')

[fdopp, fs ] = dopp_vs_theta(altitude, elevation, VG, VR, x1, y1);

%
% Figure 1- coordinate system definition
%
figure(1)

hold off
clf
plot(x, y, 'k--', x1, y1, 'k', x2, y2, 'k', xs, 0, '*k',  0, 0, 'ok', ...
     xpolar1, ypolar1, 'k', xpolar2, ypolar2, 'k',[-500 4500], [0 0], '--k', ...
     'MarkerSize',10)
 hold on
 set(gcf, 'PaperPosition', [0 0 6 4]);  
 set(gcf, 'PaperSize', [6 4])
 
 for klabel =1:size(x,2)
   ytxt = interp1(x(:,klabel), abs(y(:,klabel)), 1700) + 50;
   if( (ytxt > 0) && (delays(klabel) ~= delay1))
      text(1700+0, ytxt, ['\eta=', num2str(delaychips(klabel)),' chip'])
   end
 end
 
axis('equal')
axis([-500 4000 -2000 2000])

xlabel('x (m)')
ylabel('y (m)')
grid
hold off

%saveas(gcf, [latexfiledir, 'chgofvar'], 'pdf')

fprintf('Figure 1: \n')
fprintf('  Solid line (tau):   %2.3f chips \n', delay1chips)
fprintf('  Between lines (delta tau):   %2.3f m \n', ddelay )
fprintf('  Azimuth (theta):   %4.3f deg \n', theta1deg)
fprintf('  Delta theta:   %2.3f deg \n', dthetadeg )
fprintf('  Geometry from Example:   %2i \n', caseno )

%
% Figure 2 - Doppler illustration
%
figure(2)


hold off
clf
plot( x1, y1, 'k',  xs, 0, '*k', 0, 0, 'ok', 'MarkerSize',10)
% plot(xd, yd, 'k')
hold on

for k=1:ndopp
    [xplotd, indexd] = sort(xd(gotsoln(:,k),k));  % this makes sure that they are in order of increasing x - no overlaps
    yplotd = yd(gotsoln(:,k),k);
    plot( xplotd, yplotd(indexd), 'k--')
%    ytxt = interp(yplotd, 3500)+500;
    ytxt = interp1(xplotd, yplotd(indexd), 3500)+100;
    if ((2* floor(k/2) == k) && abs(ytxt) < 2000)
      text(3500+0, ytxt, [num2str(dopps(k)),' Hz'])
    end
 %   k
 %
end

plot( delaydopp_intersect_x, delaydopp_intersect_y, 'ok', 'MarkerSize',10, 'MarkerFaceColor',[0 0 0])

plot([-500 4500], [0 0], '--k')

ytxt = interp1(x1, abs(y1), 1700) + 100
text(1700+0, ytxt, ['\eta=', num2str(delay1chips),' chip'])
     
axis('equal')
axis([-500 4000 -2000 2000])
%    pause


xlabel('x (m)')
ylabel('y (m)')
grid
hold off


 set(gcf, 'PaperPosition', [0 0 6 4]);  
 set(gcf, 'PaperSize', [6 4])

%saveas(gcf, [latexfiledir, 'fvstheta1'], 'pdf')

fprintf('Figure 2: \n')
fprintf('  Solid line (tau):   %2.3f chips \n', delay1chips)
fprintf('  f-tilde:   %8.4f Hz \n', ftilde )
fprintf('  Geometry from Example:   %2i \n', caseno )

%
% Figure 3 - Doppler vs. theta
%
figure(3)

delaydopp_theta_intersect = atan2(delaydopp_intersect_y, delaydopp_intersect_x - x10);


plot(theta_isorange1*180/pi, fdopp-fs-ftilde, 'k', [-180 180], [0 0], 'k--', ...
    delaydopp_theta_intersect*180/pi, zeros(size(delaydopp_theta_intersect)), 'ok', ...
    'MarkerSize',10, 'MarkerFaceColor',[0 0 0])

axis([-180 180 -450 250])
xlabel('\theta (deg)')
ylabel('f_D-f (Hz)')
set(gca, 'XTick', [ -180 -135 -90 -45 0 45 90 135 180])

grid

 set(gcf, 'PaperPosition', [0 0 6 4]);  
 set(gcf, 'PaperSize', [6 4])

%saveas(gcf, [latexfiledir, 'fvstheta2'], 'pdf')

fprintf('Figure 3: \n')
fprintf('  Solid line (tau):   %2.3f chips \n', delay1chips)
fprintf('  f-tilde:   %8.4f Hz \n', ftilde )
fprintf('  Geometry from Example:   %2i \n', caseno )

%
% Fig 4 - test  - compute the Doppler along the iso-dopp lines, should be
% constant
%
%realdopp = logical(imag(xd1)==0);
%[fdopp, fs ] = dopp_vs_theta(altitude, elevation, VG, [VR(1), VR(2), 0], xd1(realdopp), yd1(realdopp));
%figure(4)

%plot(thetad1(realdopp), fdopp-fs-dopp1, '-')
