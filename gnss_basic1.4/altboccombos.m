%
% Hack to verify that the constellation in figure 2 of Rebeyrol etal 2007
% is correct for non-constant envelope boc codes.
%
% All possible combinations of the 4 spreading/data values:
%
cu = [ones(8,1); -1*ones(8,1)];
cuprime = [ones(4,1); -1*ones(4,1); ones(4,1); -1*ones(4,1)];
cl = [ones(2,1); -1*ones(2,1); ones(2,1); -1*ones(2,1); ...
      ones(2,1); -1*ones(2,1); ones(2,1); -1*ones(2,1)];
clprime = [1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1; 1; -1];
fprintf('Modulation values: \n')
cmatrix = [cu cuprime cl clprime cu+cl cuprime-clprime cuprime+clprime cu-cl]
xset = []
for s = [-1,1]
  for c = [-1, 1]
     x = ((cu + cl) * c - (cuprime - clprime)*s ) + ...
         i*((cuprime+clprime)*c + (cu-cl)*s);
     fprintf(' sine part = %2i   cosine part = %2i \n', s, c);
     x
     xset = [xset; x];
  end
end

