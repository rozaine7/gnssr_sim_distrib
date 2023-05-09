function [corrspec_fix, bandwidth_fix, bandwidth, S0] = specfix( ftilde, corrspec, Bmin)
%
% Sets a maximum correlation time on spectrum - to account for
% decorrelation due to ocean surface (breakdown of frozen-sea assumption)
%

[nlags, nfreq] = size(corrspec);

near0 = 3; % don't use S(0) due to noise, take average of 2*near0+1 points near to 0

%bandwidth = zeros(nspec,1);
%corrspec_fix = zeros(nlags,nfreq);

S0 = mean(corrspec(:,(nfreq+1)/2-near0:(nfreq+1)/2+near0),2);

bandwidth = 0.5*trapz(ftilde, corrspec, 2)./S0;
   
corrspec_fix = corrspec;

if (Bmin > 0)
    
   fprintf('Correcting spectrum to minimum bandwidth of %12.4f Hz \n', Bmin)

%for k=1:nspec
%  parest = lsqnonlin('exp2fit', [0.006, 80], [], [], [], ftilde', corrspec(k,:)');
%  bandwidth(k) = parest(2);
%  S0(k) = parest(1);
%end

   bandwidth_fix = bandwidth;
   bandwidth_fix(logical(bandwidth < Bmin)) = Bmin;

   corrspec_exp2 = S0*ones(1,nfreq).* exp(-((ones(nlags,1)*ftilde)./(bandwidth_fix*ones(1,nfreq))).^2);


   corrspec_fix(logical(bandwidth < Bmin),:) = corrspec_exp2(logical(bandwidth < Bmin),:);

else
    fprintf('No minimum bandwidth set \n') 
    bandwidth_fix = -999;
    
end

return

