function binarray = oct2bin( octarray);
%
% converts an array of octal numbers to binary.  Both input and output
% are column vectors
%

oct_table = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
oct_index = octarray+1;

binmatrix = oct_table(oct_index,:);
binarray = reshape(binmatrix',1,3*size(octarray,1))';

return
