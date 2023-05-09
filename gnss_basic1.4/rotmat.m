function T = rotmat(ax, ang)
%
% Rotation matrix about axis ax by angle ang
% angle in radians
%
switch ax
 case 1
  T = [1 0 0; 0 cos(ang) sin(ang); 0 -sin(ang) cos(ang)];
 case 2
  T = [ cos(ang) 0 -sin(ang); 0 1 0; sin(ang) 0 cos(ang)];
 case 3
   T = [ cos(ang) sin(ang) 0; -sin(ang) cos(ang) 0; 0 0 1];
 otherwise 
   display ('error axis should be 1,2, or 3')
end

return