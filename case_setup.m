%
% Loads the POLSCAT winds for comparison - use this to identify cases to
% run
%
% 11/12 by JLG for CYGNSS simulator development and TGARS paper.  Started
% with code from Justin's work in GRSL paper. 
%
% Updated with appropriate directory name for TGARS paper by JLG 1/15
%

dirname =  '/Users/jgarriso/jgarriso/projects/gnssr/airborne_simulator_paper/airborne_sim_013115/';
prefix = '2Mar09_'

eval(['load ',dirname, prefix,'level1B']);  % This is just to get elevations and prns 

select = [50 150 375 525];  % Selection of the times for comparison

pals=textread('march_3rd_to_john_with_PALS_Tb.txt');
tbv=pals(:,8);
tbh=pals(:,9);
t_pol=pals(:,1);
wpol=pals(:,4);
for j=1:length(t_pol)
    if(t_pol(j)<12)
        t_pol(j)=t_pol(j)+24;
    end
end

%t_pol=t_pol+24;

tmid_hr=t_pol*3600;

tmid = tmid_hr(select); 

%t_inds=interp1(t,[1:length(t)],t_pol,'nearest');

subplot(2,1,1)
[axh, h1, h2] = plotyy(t_pol, wpol, t_pol, tbv)


axis(axh(1), [20 26 0 40])
axis(axh(2), [20 26 120 130])
xlabel('Time (Hr)')
ylabel(axh(1),'POLSCAT Wind (m/s)')
ylabel(axh(2),'PALS   T_{b,V} (K)')
legend('POLSCAT', 'PALS')
hold on
plot(t_pol(select), wpol(select), 'o')

hold off

subplot(2,1,2)
[axh, h1, h2] = plotyy(twf/3600, prn, twf/3600, EL)
axis(axh(1), [20 26 0 40])
axis(axh(2), [20 26 0 90])
xlabel('Time (Hr)')
ylabel(axh(1),'PRN')
ylabel(axh(2),'Elevation (deg)')
legend('PRN', 'Elevation')

hold on
plot(t_pol(select), 20, 'o')

hold off

%subplot(4,1,2)
% [axh, h1, h2] = plotyy(t(1:npoints)/3600, EL*180/pi, t(1:npoints)/3600, AZ*180/pi);
%axis(axh(1), [tplotmin/3600, tplotmax/3600 0 90])
%legend('EL','AZ')
%xlabel('Time (Hr) UTC')
%ylabel(axh(1),'EL (deg)')
%ylabel(axh(2),'AZ (deg)')
%axis(axh(1), [tplotmin/3600, tplotmax/3600 0 90])
%axis(axh(2), [tplotmin/3600, tplotmax/3600 -180 180])

save tmid tmid select wpol tbv t_pol
