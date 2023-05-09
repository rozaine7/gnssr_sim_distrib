function [timetag, svn, obs1, obs2, obs3, obs4, obs5, obs6, obs7] = ...
    rinex_inpt( rinexfile, nobstypes, nobs, nheader, multiline)
%
% hack program to read in a rinex file ...
%
% multiline = 2 - will work with the MITEL C/A code only RINEX files.
%
%%rinexfile = 'dsrc3080.02o';
rid = fopen(rinexfile,'r');
%%nobstypes = 7;
%%nobs = 2400;

for i=1:nheader
  headerline = fgetl(rid);
end
svn = zeros(nobs,12);
obs1 = zeros(nobs,31);
obs2 = zeros(nobs,31);
obs3 = zeros(nobs,31);
obs4 = zeros(nobs,31);
obs5 = zeros(nobs,31);
obs6 = zeros(nobs,31);
obs7 = zeros(nobs,31);
timetag = zeros(nobs,1);
for kobs=1:nobs
  obsheader = fgetl(rid)
  YY = str2num(obsheader(1:3));
  MM = str2num(obsheader(4:6));
  DD = str2num(obsheader(7:9));
  HH = str2num(obsheader(10:12));
  MN = str2num(obsheader(13:15));
  scn = str2double(obsheader(16:26));
  eflag = str2num(obsheader(27:29));
  nsv = str2num(obsheader(30:32))
  timetag(kobs,1) = HH*3600 + MN*60 + scn;
  for ksv = 1:nsv
    svn(kobs,ksv) = str2num(obsheader(34+(ksv-1)*3:33+(ksv-1)*3+2));
    obsline1 = fgetl(rid);
    obsline1 = [obsline1, blanks(80-size(obsline1,2))];
    if (multiline == 2)
      obsline2 = fgetl(rid);
      obsline2 = [obsline2, blanks(80-size(obsline2,2))];
    end
    obs1(kobs,svn(kobs,ksv)) = str2double(obsline1(1:15));
    obs2(kobs,svn(kobs,ksv)) = str2double(obsline1(17:32));
    obs3(kobs,svn(kobs,ksv)) = str2double(obsline1(34:48));
    if (multiline == 2)
      obs4(kobs,svn(kobs,ksv)) = str2double(obsline1(50:64));
      obs5(kobs,svn(kobs,ksv)) = str2double(obsline1(65:78));
      obs6(kobs,svn(kobs,ksv)) = str2double(obsline2(1:15));  
      obs7(kobs,svn(kobs,ksv)) = str2double(obsline2(17:32));
    end
  end
end

fprintf('Read observation data from RINEX file name:  %s\n', rinexfile);

return
