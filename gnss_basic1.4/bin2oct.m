function octarray = bin2oct(binarray)
%
%  Function generates array of octal digits from array of binary digits
%  Both are column vectors
%
octal_table(:,:,1) = [0, 2; 4, 6];
octal_table(:,:,2) = [1, 3; 5, 7];
nbinary = size(binarray,1);  
npadded = 3 * ceil(nbinary/3)

binarray_padded =  [zeros(npadded-nbinary,1);binarray];

octal_sub = reshape(binarray_padded,3,npadded/3)'+1;
octal_ind = sub2ind([2,2,2], octal_sub(:,1), octal_sub(:,2), octal_sub(:,3));
octarray = octal_table( octal_ind);

return
