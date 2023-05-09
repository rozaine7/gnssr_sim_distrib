[f,Sf] = boc_spec(10,5,30e6);
deltaf = f(2)-f(1);
% Each pass band in Hz
%Filter_bands = [-2e6 -0.75e6; -0.75e6  0; ...
%                           0  0.75e6;  0.75e6 2e6]; 

taumax = 1/deltaf;
%deltatau = 1/max(f);
tau = linspace(-taumax/2, taumax/2, size(f,2));

Filter_bands = [-24e6 24e6];  % Bandlimited ACF
%Filter_bands = [-15e6 -5e6; 5e6 15e6];    % F&B sideband
%Filter_bands = [-20e6 -10e6; -10e6 0; 0 10e6; 10e6 20e6];  % 4 filters sideband

%Filter_bands = [-20e6 -15e6;  -15e6 -10e6;  ...
%                -10e6   -5e6; -5e6 0;    ...
%                 0  5e6; 5e6 10e6;  ...
%                10e6 15e6;  15e6 20e6];  
% 8 filters
 
%Filter_bands = [-20e6 -15e6;  -15e6 -10e6;  ...
%                -10e6   -5e6; -5e6 0;    ...
%                 0  5e6; 5e6 10e6;  ...
%                 15e6 20e6];   % 7 filters - representing jamming removal 10-15

Nb = size(Filter_bands,1);

Sfiltmat = [];
Acffiltmat = [];

for kb= 1:Nb
  Sfilt = Sf;
  Sfilt( logical(f < Filter_bands(kb,1))) = 0;
  Sfilt( logical(f > Filter_bands(kb,2))) = 0;
  Acffilt = ifftshift(ifft(fftshift(Sfilt)));
  Sfiltmat = [Sfiltmat; Sfilt];
  Acffiltmat = [Acffiltmat; Acffilt];
end
  

