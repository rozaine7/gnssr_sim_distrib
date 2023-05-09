function eph = yuma_read( ephemfile)
%
% Read a text Yuma-format ephemeris file.
%
fid = fopen(ephemfile, 'r');
eph = [];
eofflag =0;

while (eofflag == 0)
  ephsv = zeros(23,1);
  firstline = fgetl(fid);
  if ( size(firstline,2) > 4)
    ephline = fgetl(fid)
    ephsv(1,1) = str2num(ephline(29:end));
    fprintf(' Reading Ephemeris for PRN %4i \n', ephsv(1,1));
    for k=2:17
      ephline = fgetl(fid);
      ephsv(k,1) = str2num( ephline(29:45));
    end
    ephline = fgetl(fid);
    ephsv(18,1) = str2num(ephline(29:end));
    for k=19:23
      ephline = fgetl(fid);
      ephsv(k,1) = str2num( ephline(29:45));
    end
    skipline = fgetl(fid);
    eph = [eph, ephsv];
  else
    eofflag = 1;
  end
end



return