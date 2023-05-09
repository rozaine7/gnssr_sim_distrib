function [re] = eci2ecef(ri,t,theta0, we) ; 
% convert eci coordinates to ecef coordinates
%global we - don't use global variables in this version -JLG
%
if nargin < 4
    sidday = 23*3600+56*60+4.090524 ; % sec
    we = 2*pi / sidday ;  % 1/sec
end

theta = theta0 + we * t ; 

re = zeros(size(ri)) ; 

% T * r.'
re(:,1) =  cos(theta) .* ri(:,1) + sin(theta) .* ri(:,2) ; 
re(:,2) = -sin(theta) .* ri(:,1) + cos(theta) .* ri(:,2) ; 
re(:,3) =  ri(:,3) ; 

% T * [v - wxr].' 
if size(ri,2) == 6 | size(ri,2) == 9 
   re(:,4) =  cos(theta) .* (ri(:,4) + we*ri(:,2)) ...
            + sin(theta) .* (ri(:,5) - we*ri(:,1)) ; 
   re(:,5) = -sin(theta) .* (ri(:,4) + we*ri(:,2)) ...
            + cos(theta) .* (ri(:,5) - we*ri(:,1)) ; 
   re(:,6) =  ri(:,6) ; 
end

% T * [a - wxwxr - 2*wxr].'
if size(ri,2) == 9 
   re(:,7) = cos(theta).*(ri(:,7) + we^2*ri(:,1) + 2*we*ri(:,5)) ...
           + sin(theta).*(ri(:,8) + we^2*ri(:,2) - 2*we*ri(:,4)) ; 
   re(:,8) =-sin(theta).*(ri(:,7) + we^2*ri(:,1) + 2*we*ri(:,5)) ...
           + cos(theta).*(ri(:,8) + we^2*ri(:,2) - 2*we*ri(:,4)) ; 
   re(:,9) =  ri(:,9) ; 
end
