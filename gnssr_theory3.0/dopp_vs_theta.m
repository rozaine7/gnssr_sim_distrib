function [fdopp, fs] = ...
         dopp_vs_theta( alt, gammadeg, VGPS, VR, xisor, yisor)
%
% Corrected 2/15 by JLG.
% compute Doppler around a given iso-range ellipse
% Assuming a flat Earth - trasnmitter at infinity
%
% df is the set of Doppler frequnecies, in Hz, expected to be in a row
% vector. Columns of outputs all correspond to a single iso-Dopp.
%
% Polar coordinates are used in this calculation, except the origin in at
% the x=0, not at the specular point
% 
% Run isorange first - to generate the xisor and yisor coordinates on the
% surface.  These must be column vectors for a single iso-range ellipse. 
% 
gamma = gammadeg * pi/180;
lambda = 0.19;

ns = [-cos(gamma); 0; sin(gamma)];  % unit vector in direction of specular reflection
m = [-cos(gamma); 0; -sin(gamma)]; % unit vector in direction of incident radiation
fs = -(dot(ns,VR) - dot(m,VGPS))/lambda; % Doppler at specular point

%size(xisor)
%size([xisor  yisor 0*ones(size(xisor,1),1)])
%size([-alt*cos(gammadeg*180/pi) 0 alt*sin(gammadeg*180/pi)])

Rvec = ones(size(xisor,1),1)*[0 0 alt] - [xisor  yisor 0*ones(size(xisor,1),1)];
Rmag = sqrt( Rvec(:,1).^2 + Rvec(:,2).^2 + Rvec(:,3).^2);

n = Rvec ./ (Rmag*ones(1,3));

size(n)
%size( dot(n,VR,2))

fdopp = -(dot(n,(ones(size(xisor,1),1)*VR),2) - dot(m,VGPS))/lambda;  % doppler at each point (xisor, yisor)



return