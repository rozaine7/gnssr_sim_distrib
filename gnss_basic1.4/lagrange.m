function xt = lagrange(t, x0, t0)
%
% Lagrange interpolation at time t, from series of data x0 at time t0
%
[ndata, xdim] = size(x0);
xt = 0;

for i=1:ndata
  ln = 1;
  ld = 1;
  for k=1:ndata
    if (k ~= i)
       ln = ln*(t - t0(k));
       ld = ld*(t0(i) - t0(k));
    end
  end
  xt = xt + x0(i,:)*(ln/ld);
end


return
