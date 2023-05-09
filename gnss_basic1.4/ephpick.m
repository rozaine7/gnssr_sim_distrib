function ephvec = ephpick(ephmat, prn, gpstime);
%
% selects the ephemeris with the correct prn, and the TOE cloeset to 
% gpstime - ephmat should be in the Yuma format, with 1 colm per sat/per
% time
%
correctsv = logical( ephmat(1,:) == prn);
ephsv = ephmat(:,correctsv);

[mintime, ptr ] = min( abs(ephsv(3,:) - gpstime));

ephvec = ephsv(:,ptr);

return