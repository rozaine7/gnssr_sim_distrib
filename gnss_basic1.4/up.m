function p = up(x)
%
% Unit pulse function
%
p = zeros(size(x));
p(logical(abs(x)<1)) = 1;
return;