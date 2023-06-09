function eph = rinexe(ephemerisfile)
%
% Reads a RINEX Navigation Message file from Builder-2
%   	  reformats the data into a matrix with 21
%	     rows and a column for each satellite.
%	     The matrix is stored in outputfile
%

% Typical call: rinexe('pta.96n','pta.nav')

%
% Based upon original code:
%Kai Borre 04-18-96
%Copyright (c) by Kai Borre
%$Revision: 1.0 $  $Date: 1997/09/24 $

%Units are either seconds, meters, or radians
fide = fopen(ephemerisfile);
head_lines = 0;
while 1			    % We skip header
   head_lines = head_lines+1;
   epline = fgetl(fide);
   answer = findstr(epline,'END OF HEADER');
   if  ~isempty(answer), break;  end;
end;
head_lines
noeph = -1;
while 1
   noeph = noeph+1;
   epline = fgetl(fide);
   if epline == -1, break;  end
end;
noeph = noeph/8
frewind(fide);
for i = 1:head_lines, epline = fgetl(fide); end;

% Set aside memory for the input
svprn = zeros(1,noeph);
weekno = zeros(1,noeph);
toc   = zeros(1,noeph);
tgd   = zeros(1,noeph);
aodc  = zeros(1,noeph);
toe   = zeros(1,noeph);
af2   = zeros(1,noeph);
af1   = zeros(1,noeph);
af0   = zeros(1,noeph);
aode  = zeros(1,noeph);
deltan = zeros(1,noeph);
M0    = zeros(1,noeph);
ecc   = zeros(1,noeph);
roota = zeros(1,noeph);
toe   = zeros(1,noeph);
cic   = zeros(1,noeph);
crc   = zeros(1,noeph);
cis   = zeros(1,noeph);
crs   = zeros(1,noeph);
cuc   = zeros(1,noeph);
cus   = zeros(1,noeph);
Omega0 = zeros(1,noeph);
omega = zeros(1,noeph);
i0    = zeros(1,noeph);
Omegadot = zeros(1,noeph);
idot  = zeros(1,noeph);
accuracy = zeros(1,noeph);
health = zeros(1,noeph);
fit = zeros(1,noeph);

for i = 1:noeph
   epline = fgetl(fide);	  %
   size(epline)
   svprn(i) = str2num(epline(1:2));
   year = epline(3:6);
   month = epline(7:9);
   day = str2num(epline(10:12));
   hour = str2num(epline(13:15));
   minute = str2num(epline(16:18));
   second = str2num(epline(19:22));
%   toc(i) = day*3600*24 + hour*3600 + minute*60 + second;
   af0(i) = str2num(epline(23:41));
   af1(i) = str2num(epline(42:60));
   af2(i) = str2num(epline(61:79));
   epline = fgetl(fide);	  %
   IODE = epline(5:22);
   crs(i) = str2num(epline(23:41));
   deltan(i) = str2num(epline(42:60));
   M0(i) = str2num(epline(61:79));
   epline = fgetl(fide);	  %
   cuc(i) = str2num(epline(5:22));
   ecc(i) = str2num(epline(23:41));
   cus(i) = str2num(epline(42:60));
   roota(i) = str2num(epline(61:79));
   epline=fgetl(fide);
   toe(i) = str2num(epline(5:22));
   toc(i) = toe(i);
   cic(i) = str2num(epline(23:41));
   Omega0(i) = str2num(epline(42:60));
   cis(i) = str2num(epline(61:79));
   epline = fgetl(fide);	    %
   i0(i) =  str2num(epline(5:22));
   crc(i) = str2num(epline(23:41));
   omega(i) = str2num(epline(42:60));
   Omegadot(i) = str2num(epline(61:79));
   epline = fgetl(fide);	    %
   idot(i) = str2num(epline(5:22));
   codes = str2num(epline(23:41));
   weekno = str2num(epline(42:60));
   L2flag = str2num(epline(61:79));
   epline = fgetl(fide);	    %
   svaccur = str2num(epline(5:22));
   svhealth = str2num(epline(23:41));
   tgd(i) = str2num(epline(42:60));
   iodc = epline(61:79);
   epline = fgetl(fide);	    %
   tom(i) = str2num(epline(5:22));
%   spare = epline(23:41);
%   spare = epline(42:60);
%   spare = epline(61:79);
end
status = fclose(fide)

%  Description of variable eph. 
eph(1,:) = svprn;
eph(2,:) = af2;
eph(3,:) = M0;
eph(4,:) = roota;
eph(5,:) = deltan;
eph(6,:) = ecc;
eph(7,:) = omega;
eph(8,:) = cuc;
eph(9,:) = cus;
eph(10,:) = crc;
eph(11,:) = crs;
eph(12,:) = i0;
eph(13,:) = idot;
eph(14,:) = cic;
eph(15,:) = cis;
eph(16,:) = Omega0;
eph(17,:) = Omegadot;
eph(18,:) = toe;
eph(19,:) = af0;
eph(20,:) = af1;
eph(21,:) = af2;
eph(22,:) = toc;

%fidu = fopen(outputfile,'w');
%count = fwrite(fidu,[eph],'double');
fclose all
%%%%%%%%% end rinexe.m %%%%%%%%%
