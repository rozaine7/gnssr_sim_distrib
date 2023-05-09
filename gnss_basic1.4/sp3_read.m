function [gpsxyz, cbias, gpst] = sp3_read( sp3file )
%
% Read a text ephemeris data in sp3_read.
%
sp3id = fopen(sp3file);

%for i=1:22
%  headerlines = fgetl(sp3id);
%end

headerline1 = fgetl(sp3id);
headerline2 = fgetl(sp3id);
headerline3 = fgetl(sp3id);

for i=4:22
  headerlines = fgetl(sp3id);
end

nsv = str2num(headerline3(5:6));

ntimes = str2num(headerline1(38:39));

svx = zeros(ntimes,31);
svy = zeros(ntimes,31);
svz = zeros(ntimes,31);
cbias = zeros(ntimes,31);
svt = zeros(ntimes,1);

for ktime = 1:ntimes
  dateline = fgetl(sp3id);
  YY = str2num(dateline(4:7));
  MM = str2num(dateline(9:10));
  DD = str2num(dateline(12:13));
  HH = str2num(dateline(15:16));
  MN = str2num(dateline(18:19));
  scn = str2double(dateline(22:31));
  svt(ktime,1) = HH*3600 + MN*60 + scn;
  for ksv = 1:nsv
    satline = fgetl(sp3id);
    if (satline(1) == 'P')
      svid = str2num(satline(3:4));
      svx(ktime,svid) = str2double(satline(6:18));
      svy(ktime,svid) = str2double(satline(20:32));
      svz(ktime,svid) = str2double(satline(34:46));
      cbias(ktime,svid) = str2double(satline(48:60));
    end
  end
end

gpst = svt;
gpsxyz(:,1,:) = svx * 1000;
gpsxyz(:,2,:) = svy * 1000;
gpsxyz(:,3,:) = svz * 1000;

fprintf('Read from Precise Ephemeris File: %s \n', sp3file);
return