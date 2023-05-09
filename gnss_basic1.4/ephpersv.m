function [ephsort, ageofdata] = ephpersv( ephemblock, midtime)
%
% Extracts ephemeris data from the format provided by the Misra and Enge
% books.  Pulls out ephemeris closest to midtime
%
prnset = ephemblock(:,2);
eph = zeros(24,32);
ephsort = zeros(24,32);
for k=1:32
  ephset = ephemblock(logical(prnset == k),:);
  if (size(ephset,1) > 0)
    [val, ptr] = min(abs(midtime - ephset(:,1)));
    eph(:,k) = ephset(ptr,:)';
    ageofdata(1,k) = val;
  end
end

ephsort(1,:) = eph(2,:);
ephsort(2,:) = eph(9,:);
ephsort(3,:) = eph(4,:);
ephsort(4,:) = eph(15,:);
ephsort(5,:) = eph(16,:);
ephsort(6,:) = eph(10,:);
ephsort(7,:) = eph(14,:);
ephsort(8,:) = eph(13,:);
ephsort(9,:) = eph(12,:);
ephsort(10,:) = eph(11,:);
ephsort(11,:) = eph(17,:);
ephsort(12,:) = eph(19,:);
ephsort(13,:) = eph(18,:);
ephsort(14,:) = eph(23,:);
ephsort(15,:) = eph(22,:);
ephsort(16,:) = eph(21,:);
ephsort(17,:) = eph(20,:);
ephsort(20,:) = eph(3,:);
ephsort(21,:) = eph(5,:);
ephsort(22,:) = eph(6,:);
ephsort(23,:) = eph(7,:);
ephsort(24,:) = eph(8,:);

return

