function S2 = elfouhaily2d (k, phi, U10, inverse_wave_age)
%
% 2d wave spectrum from Elfouhaily - takes a vector k and scalar phi
%
S1 = elfouhaily1d (k, U10, inverse_wave_age);
cm = 0.23;
g = 9.81;
km = 370;
ustar = 0.023 *  U10^(1.23);
a0 = log(2)/4;
ap = 4;
am = 0.13 * ustar/cm;
cp = U10 / inverse_wave_age;
c = sqrt(g * (1+(k/km).^2) ./ k);
Delta = tanh( a0 + ap * ( c/cp).^(2.5) + am * (cm ./ c).^(2.5));
spreading_function = (1/(2*pi)) * ( 1 + Delta .* cos(2 * phi));
S2 = S1 .* spreading_function ./ k;

return