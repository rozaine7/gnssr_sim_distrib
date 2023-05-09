function doy = YMD2DOY(Y, M, D);
%
% DOY = MD2DOY( Y, M, D)
% 
% Computes the day of year (DOY) from the year (Y), month (M) and day (D).  
% Y, M and D can be matrices (same dimensions)
%

  doy = datenum(Y, M, D) - datenum(Y, 1,1)+1;

return
