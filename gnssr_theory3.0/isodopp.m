function [x, y, theta, beta, fs, gotsoln,r ] = ...
         isodopp( alt, gammadeg, df, VGPS, VR, maxtheta, lambda, maxwidth)
%
% Mod 5/2016 y JLG to remove hard-coding at GPS L1 - fixed (hack) the
% problem with zenith geometry - maxwidth is only used for that special
% case.
%
% Corrected 10/14 by JLG.
% compute bistatic iso-Doppler lines assuming a flat Earth - trasnmitter at
% infinity
%
% df is the set of Doppler frequnecies, in Hz, expected to be in a row
% vector. Columns of outputs all correspond to a single iso-Dopp.
%
% Polar coordinates are used in this calculation, except the origin in at
% the x=0, not at the specular point
% 
gamma = gammadeg * pi/180;

ndopp = size(df,2);



if (nargin < 7)
  lambda = 0.19;
end

ns = [-cos(gamma); 0; sin(gamma)];  % unit vector in direction of specular reflection
m = [-cos(gamma); 0; -sin(gamma)]; % unit vector in direction of incident radiation
fs = -(dot(ns,VR) - dot(m,VGPS))/lambda; % Doppler at specular point
fdopp = fs + df;
C2 = lambda * fdopp - dot(m,VGPS);

theta = -maxtheta:0.001:maxtheta;  % increased step size from 0.0001 2/2015 jlg
theta = theta';


ntheta = size(theta,1)


%
% revised, simpler version - will use principle that \hat{n} dot V^R = -
% C2
cosbeta = (ones(ntheta,1)*C2) ./ ((cos(theta)*VR(1) + sin(theta)*VR(2))*ones(1,ndopp));
beta = acos(cosbeta);
gotsoln = logical(imag(beta) == 0);


r = alt./tan(beta);
x = r.*cos(theta*ones(1,ndopp));
y = r.*sin(theta*ones(1,ndopp));


if (gammadeg == 90)
    fprintf('WARNING- specular iso-Dopp is not generated for zenith (gamma = 90)\n')
    middop = (ndopp+1)/2;
    if (VR(2) ~= 0)
       x(:,middop) = linspace(-maxwidth, maxwidth, ntheta)';
       y(:,middop) = -VR(1)/VR(2)* x(:,middop);
    else
       y(:,middop) = linspace(-maxwidth, maxwidth, ntheta)';
       x(:,middop) = -VR(2)/VR(1)* x(:,middop);
    end
end

return