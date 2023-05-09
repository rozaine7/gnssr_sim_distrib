function out = constellation_plot(s, sym)
%
% Generate a constellation plot for a complex modulation
%
npts = size(s)
hold on
for k=1:npts
  plot(s(k), sym)
  text(real(s(k))+0.1, imag(s(k))+0.1, num2str(k))
end