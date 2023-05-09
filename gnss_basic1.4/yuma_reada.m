function alm = yuma_read( almfile)
%
% Read a text Yuma-format ephemeris file.
%
fid = fopen(almfile, 'r');
alm = [];
eofflag =0;

    % tabs got all messed up - will read each line separately, only dependent
    % upon the length of the title string.
    element_name_length(1) = size(char('ID:'),2);
    element_name_length(2) = size(char('Health:'),2);
    element_name_length(3) = size(char('Eccentricity:'),2);
    element_name_length(4) = size(char('Time of Applicability(s):'),2);
    element_name_length(5) = size(char('Orbital Inclination(rad):'),2);
    element_name_length(6) = size(char('Rate of Right Ascen(r/s):'),2);
    element_name_length(7) = size(char('SQRT(A)  (m 1/2):'),2);
    element_name_length(8) = size(char('Right Ascen at Week(rad):'),2);
    element_name_length(9) = size(char('Argument of Perigee(rad):'),2);
    element_name_length(10) = size(char('Mean Anom(rad):'),2);
    element_name_length(11) = size(char('Af0(s):'),2);
    element_name_length(12) = size(char('Af1(s/s):'),2);
    element_name_length(13) = size(char('week:'),2);

while (eofflag == 0)
  almsv = zeros(13,1);
  firstline = fgetl(fid);
  if ( size(firstline,2) > 4)

    for k=1:13
       almline = fgetl(fid);
       almsv(k,1) = str2num(almline(element_name_length(k)+1:end));
    end
    fprintf(' Read Almanac for PRN %4i week %5i \n', almsv(1,1), almsv(13,1));
    skipline = fgetl(fid);
    alm = [alm, almsv];
  else
    eofflag = 1;
  end
end



return