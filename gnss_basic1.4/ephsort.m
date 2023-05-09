function esortmat = ephsort(eph);
%
% Re-arranges elements from the data files distributed on Misra/Enge CD
% ("Format A") to the format of the Yuma ephemeris.  
%
esortmat(1,:) = eph(2,:);  % PRN
esortmat(2,:) = eph(9,:);  % eccentricity
esortmat(3,:) = eph(4,:);  % TOE
esortmat(4,:) = eph(15,:); % i0
esortmat(5,:) = eph(16,:); % rate of RAN (Omega dot)
esortmat(6,:) = eph(10,:); % sqrt(a)
esortmat(7,:) = eph(14,:); % LAN0
esortmat(8,:) = eph(13,:); % Argp
esortmat(9,:) = eph(12,:); % M0
esortmat(10,:) = eph(11,:);% n dot 
esortmat(11,:) = eph(17,:);% i dot
esortmat(12,:) = eph(19,:);% cuc
esortmat(13,:) = eph(18,:);% cus
esortmat(14,:) = eph(23,:);% crc
esortmat(15,:) = eph(22,:);% crs
esortmat(16,:) = eph(21,:);% cic
esortmat(17,:) = eph(20,:);% cis
esortmat(20,:) = eph(3,:); % TOC
esortmat(21,:) = eph(5,:); % af0
esortmat(22,:) = eph(6,:); % af1
esortmat(23,:) = eph(7,:); % af2

return