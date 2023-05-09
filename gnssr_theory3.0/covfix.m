function Rnonzero = covfix( R, band, threshold, evthreshold, fixflag )
%
% Removes some numerical problems in bin-bin covariance matrix. 
%
% fixflag = 1: band and thrshold only
%         = 2: normalize - make a correlation coefficient
%
ndims = size(R,1);

%
% Zero off-diagonals outside of +/- band
%
Rband = zeros(ndims,ndims);
for k=1:ndims
    Rband(k,max(1,k-band):min(ndims,k+band)) = R(k,max(1,k-band):min(ndims,k+band));
end

%
% Set minimum threshold on diagonal - this might be redundant with the
% minimum eigenvalues - but keep it in for now.
%
for k=1:ndims
    Rband(k,k) = max(threshold, Rband(k,k));
end

if fixflag == 2
% 
% Normalize to make diagonals unity 
%
    Rnorm1 = zeros(ndims,ndims);
    Rnorm2 = zeros(ndims,ndims);

    for k=1:ndims
        Rnorm1(k,:) = Rband(k,:)/sqrt(Rband(k,k));
    end

    for k=1:ndims
        Rnorm2(:,k) = Rnorm1(:,k)/sqrt(Rband(k,k));
    end

    Rfix = Rnorm2;

else
    Rfix = Rband;
end
   
%
% Fix numerical problems and force it to be symmetric
%

Rfix = (Rfix + Rfix')/2;

%
% Set minimum thershold on eigenvalues to allow Cholesky Factorization
%
[V, D] = eig(Rfix);
Dfix = zeros(ndims,1);
for k =1:ndims
    Dfix(k) = max(evthreshold, D(k,k));
end

Rnonzero = V * diag(Dfix)* V'; 


end

