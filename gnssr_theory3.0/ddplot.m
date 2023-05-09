%
% quick script to generate a isorange-isoDoppler plot
%
VR = 150*[-sin(30*pi/180) -cos(30*pi/180) 0];
VGPS = 487*[-0.53, -0.74, 0.42];
gammadeg = 51.0;
alt = 4000;
df = [ -500 -400 -300 -200 -100 0.1 100 200 300 400 500 600 700 800 800 1000];
dtau = [0.5 1 1.5 2 2.5 3 3.5 4]*293;

clf
hold

for k=1:size(df,2)
  [x, y, thet, bet] = isodopp(alt, gammadeg, df(k), VGPS, VR, pi/2);
  xpos = logical(bet > 0);
  xneg = logical(bet <= 0);
  s1 = plot(x(xpos), y(xpos), 'y');
  s2 = plot(x(xneg), y(xneg), 'y');
  set(s1, 'LineWidth', 2);
  set(s2, 'LineWidth', 2);
end


for k=1:size(dtau,2)
  [x, y, r, theta0, x0, xs] = isorange(alt, gammadeg, dtau(k));
  s3 = plot(x,y, 'y')
  set(s3, 'LineWidth', 2);
end


