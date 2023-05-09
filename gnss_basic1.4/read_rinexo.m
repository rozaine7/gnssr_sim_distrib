function [T, L1, L2, C1, P1, P2] = read_rinexo(rinexo);

% READ_RINEXO	Read RINEX observation file
%	[T, L1, L2, C1, P1, P2] = read_rinexo(rinexo)
%	T = time vector (hours)
%	L1, L2, C1, P1, P2 are matrices with one satellite per
%	column and epochs in the row direction.
%	Ex. L1(100,13) will give L1 phase at epoch 100 for PRN 13

%
% Modified 8/14 for AAE575 to get updated CORS files. 
%
% hard-wired for now but should be read from rinex header
nobs = 5;
samplint = 30;
total_epoch = 2880;
total_nsat = 31;

% declare observable matrices
L1=zeros(total_epoch,total_nsat);
L2=zeros(total_epoch,total_nsat);
C1=zeros(total_epoch,total_nsat);
P1=zeros(total_epoch,total_nsat);
P2=zeros(total_epoch,total_nsat);
T=zeros(total_epoch,1);

fid=fopen(rinexo,'r');
j = 0;
while 1
  line = fgetl(fid);
  if ~isstr(line), break, end
  if (length(line)>32 & line(33)=='G')
     j = j+1;
     nsat = str2num(line(31:32));
     [t] = sscanf(line, '%d %d %d %d %d %d');
     T(j) = t(4)+t(5)/60+t(6)/3600;
     sv = [];
     i = 34;
     while (i < length(line))
         sv = [sv;str2num(line(i:i+1))];
         i = i+3;
     end
     for i = 1:nsat
         line = fgetl(fid);
         if ((length(line))>=14 & ~isempty(str2num(line(1:14)))), L1(j,sv(i)) = str2num(line(1:14));, end
         if ((length(line))>=30 & ~isempty(str2num(line(17:30)))), L2(j,sv(i)) = str2num(line(17:30));, end
         if ((length(line))>=46 & ~isempty(str2num(line(33:46)))), C1(j,sv(i)) = str2num(line(33:46));, end
         if ((length(line))>=62 & ~isempty(str2num(line(49:62)))), P1(j,sv(i)) = str2num(line(49:62));, end
         if ((length(line))>=78 & ~isempty(str2num(line(65:78)))), P2(j,sv(i)) = str2num(line(65:78));, end
     end
  end
end
fclose(fid);
