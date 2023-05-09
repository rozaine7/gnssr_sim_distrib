function G = gmatrix ( rcvrxyz, gpsxyz)
%
% linearized observation equation 
%
nsats = size(gpsxyz,1);
rcvrmat = ones(nsats,1) * rcvrxyz(1:3);
los = gpsxyz - rcvrmat;
range = sqrt(sum(los.^2,2));

G = [-los ./ (range * ones(1,3)), ones(nsats,1)];

return