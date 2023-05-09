    nlagspre = size(pwr,1)
    nacf = size(acfamp,1)
    nlags = nlagspre + nacf;
    acfrow = [acfamp', zeros(1,nlagspre)];
    for k=1:nlags

       pwr3 = [zeros(1,k-1), acfrow(1:nlags-k+1)] .* [zeros(1,nacf), pwr'];

       Rtau(k,:) = delaystep * conv(acfamp, pwr3);
       pwr3mat(k,:) = pwr3;
       fprintf('%4i ', k)
    end
   Rtau = Rtau(2:end,nacf:end-1);