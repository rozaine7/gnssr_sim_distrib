function [xecef, yecef, zecef, sv_clock, trel] = ephem(eph, gpstime)
%
%  Program to read in the GPS satellite ephemeris data and
%  compute a cartesian position from it.
%
%  Converted from FORTRAN code, written for Individual GPS Lab project -
%  December, 1993 by James L. Garrison 
%
%  Input:  eph = 21 X Nsat = broadcast ephemeris for N satellites.
%          gpstime = Ntime X 1 - time desired for ephemeris is GPS time
%          sec.
%  Output: xecef, yecef, zecef, sv_clock - ea. Ntime X Nsat - ECEF
%          positons and clock offsets of GPS sats at times in gpstime.
%          sv_clock and trel are in seconds
%
%  8-22-03 = Clock bias terms included and relativity correction
%  included.  Allowed matrix for gpstime to be used - allowing a
%  different time for computation for each GPS satellite.
%
mu = 3.986005e14;
omega_earth = 7.2921151467e-5;

%
% Defined Ephemeris data are read from input matrix
%
nsats    = size(eph,2);
ntimes   = size(gpstime,1);
prn      = eph(1,:);
mean_0   = eph(9,:);
sqrta    = eph(6,:);
deltan   = eph(10,:);
eccen    = eph(2,:);
argp     = eph(8,:);
cuc      = ones(ntimes,1)*eph(12,:);
cus      = ones(ntimes,1)*eph(13,:);
crc      = ones(ntimes,1)*eph(14,:);
crs      = ones(ntimes,1)*eph(15,:);
inc_0    = ones(ntimes,1)*eph(4,:);
inc_dot  = ones(ntimes,1)*eph(11,:);
cic      = ones(ntimes,1)*eph(16,:);
cis      = ones(ntimes,1)*eph(17,:);
node_0   = ones(ntimes,1)*eph(7,:);
node_dot = ones(ntimes,1)*eph(5,:);
t_0      = eph(3,:);
t_c      = ones(ntimes,1)*eph(20,:);
af0      = ones(ntimes,1)*eph(21,:);
af1      = ones(ntimes,1)*eph(22,:);
af2      = ones(ntimes,1)*eph(23,:);


%
%  Preliminary data are computed
%
      sma = sqrta.^2;
      n_0 = sqrt (mu*ones(1,nsats)./(sma.^3));      
%
%  Kepler's equation is solved by FSOLVE to obtain the
%  eccentric and then the true anomaly
%
         t_k = gpstime*ones(1,nsats) - ones(ntimes,1)*t_0;
         c_n = n_0 + deltan;
         mean_k = ones(ntimes,1)*mean_0 + (ones(ntimes,1)*c_n) .* t_k;
	 eanom_k = zeros(ntimes,nsats);
	 options = optimset('Jacobian','on','Diagnostics','off', ...
			    'Display','off','TolFun',1e-14,'TolX',1e-14);
         fprintf('  Solving Keplers equation for PRN: \n');
	 for i=1:nsats
	   if (sqrta(i) > 0)
             for k=1:ntimes
             % 
             % Simpler solution of KE - don't need that high of accuracy
             E = mean_k(k,i);
             while (abs(E-mean_k(k,i))>1e-8)
               E = mean_k(k,i) - eccen(i)*mean_k(k,i);
             end
             eanom_k(k,i) = E;
%	     eanom_k(k,i) = fsolve('kepler_eqn', mean_k(k,i), options, ...
%				   mean_k(k,i), eccen(i));
             end
	   end
           fprintf('%2i, ', i);
	 end  
         fprintf(' ... done\n');
         truan_over2 = ...
	     atan(ones(ntimes,1)*(sqrt((1 + eccen)./(1 - eccen))) ...
	      .*tan(eanom_k/2));
         truan_k = 2*truan_over2;
         r = (ones(ntimes,1)*sma) .* ...
	     (1 - (ones(ntimes,1)*eccen) .* cos(eanom_k));
%
%  The corrections are applied based upon the argument of latitude
%   
         arglat_k = (ones(ntimes,1)*argp) + truan_k;
         del_arglat = cus .* sin(2* arglat_k) + cuc .* cos(2 * arglat_k);
         del_r = crs .* sin(2 * arglat_k) + crc .* cos(2 * arglat_k);
         del_inc = cis .* sin(2 * arglat_k) + cic .* cos(2 * arglat_k); 
         c_arglat_k = arglat_k + del_arglat;
         c_r_k = r + del_r;
         c_inc_k = inc_0 + del_inc + inc_dot .* t_k;
%
%  The position of the satellite in the orbit plane is found
%
         x = c_r_k .* cos (c_arglat_k);
         y = c_r_k .* sin (c_arglat_k);
%
%  These are rotated about the line of nodes to the primed
%  coordinate system 
%

         xprime = x;
         yprime = y .* cos(c_inc_k);
         zprime = y .* sin(c_inc_k);
%
%  The prime coordinate system is rotated in the equatorial plane
%  through the longitude of the ascending node so as to put the new
%  x-axis through the Greenwich meridian
%
         long_node_k = node_0 + ...
	     (node_dot - (omega_earth*ones(ntimes,nsats))) .* t_k  ...
	     - (omega_earth*ones(ntimes,nsats)) .* (ones(ntimes,1)*t_0);
         xecef = xprime .* cos(long_node_k) -  yprime .* sin(long_node_k);
         yecef = xprime .* sin (long_node_k) +  yprime .* cos(long_node_k);
         zecef = zprime;
	 
%
%  The satellite clock correction term is computed.
%
         sv_clock = af0 + af1.*(gpstime*ones(1,nsats)-t_c) + ...
	     af2.*(gpstime*ones(1,nsats)-t_c).^2;
%
% Relativistic correction
%
         F = -4.442807633E-10;  %sec/sqrt(meter)
	 trel = F * (ones(ntimes,1)*(eccen .* sqrta)) .* sin(eanom_k);
     
	 

      return

