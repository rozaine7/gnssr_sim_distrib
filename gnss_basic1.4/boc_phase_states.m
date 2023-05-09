%
%  Generate table 8 in the Galileo ICD - for the 8-PSK phase states for
%  BOC codes.
%
%  For simplicity, frequency will be in units of fs (one subcarrier cycle
%  will be exactly 1).
%
t = [0:0.001:8]'/8;    
tprime = t - 1/4;
scs = sqrt(2)/4 * sign( cos( 2 * pi * t - pi/4)) + ...
    0.5 * sign(cos(2 * pi * t)) + sqrt(2)/4 * sign( cos( 2 * pi * t+pi/4));
scp = -sqrt(2)/4 * sign( cos( 2 * pi * t - pi/4)) +  ...
    0.5 * sign(cos(2 * pi * t)) - sqrt(2)/4 * sign( cos(2*pi*t+pi/4));

scsprime = sqrt(2)/4 * sign( cos( 2 * pi * tprime - pi/4)) + ...
    0.5 * sign(cos(2 * pi * tprime)) + sqrt(2)/4 * sign( cos( 2 * pi * tprime+pi/4));
scpprime = -sqrt(2)/4 * sign( cos( 2 * pi * tprime - pi/4)) +  ...
    0.5 * sign(cos(2 * pi * tprime)) - sqrt(2)/4 * sign( cos(2*pi*tprime+pi/4));

e5aI = [-1* ones(1,8)  ones(1,8)];
e5bI = [-1*ones(1,4)  ones(1,4) -1*ones(1,4)  ones(1,4)];
e5aQ = [-1*ones(1,2) ones(1,2) -1*ones(1,2) ones(1,2) ...
        -1*ones(1,2) ones(1,2)  -1*ones(1,2) ones(1,2) ];
e5bQ = [-1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 ];

e5aI = ones(8,1) * e5aI;
e5aQ = ones(8,1) * e5aQ;
e5bI = ones(8,1) * e5bI;
e5bQ = ones(8,1) * e5bQ;

ebar5aI = e5aQ .* e5bI .* e5bQ;
ebar5aQ = e5aI .* e5bI .* e5bQ;
ebar5bI = e5bQ .* e5aI .* e5aQ;
ebar5bQ = e5bI .* e5aI .* e5aQ;

scstab = 0.5*[sqrt(2)+1; 1; -1; -sqrt(2)-1; -sqrt(2)-1; -1; 1; sqrt(2)+1];
scptab = 0.5*[-sqrt(2)+1; 1; -1; sqrt(2)-1; sqrt(2)-1; -1; 1; -sqrt(2)+1];

index4 = mod([0:7]-2,8)+1;  % index of coefficient 2 steps back 

subcarrsc1 = scstab - i * scstab(index4');
subcarrsc2 = scstab + i * scstab(index4');

subcarrpc1 = scptab - i * scptab(index4');
subcarrpc2 = scptab + i * scptab(index4');


subcarrsc1 = subcarrsc1 * ones(1,16);
subcarrsc2 = subcarrsc2 * ones(1,16);
subcarrpc1 = subcarrpc1 * ones(1,16);
subcarrpc2 = subcarrpc2 * ones(1,16);

e5 = 1/(2 * sqrt(2)) * (   (e5aI + i*e5aQ) .* subcarrsc1 + ...
                           (e5bI + i*e5bQ) .* subcarrsc2 + ...
                      (ebar5aI + i*ebar5aQ) .* subcarrpc1 + ...
                      (ebar5bI + i*ebar5bQ) .* subcarrpc2);

k = mod(angle(e5)/(pi/4),8);  % This should give the result i table 8 of
                              % the Galileo OS ICD, with the exception
                              % that point '8' in that table will be '0'
                              % here (2pi = 0).

k(logical(k==0))=8;

%
% Write these out a latex table
%
fid = fopen('altboc_phase_states.tex', 'w');

fprintf(fid,'$e_{aI}$')
fprintf(fid, '& %3i ', e5aI(1,:));
fprintf(fid, '\\\\ \n')
fprintf(fid,'$e_{bI}$')
fprintf(fid, '& %3i ', e5bI(1,:));
fprintf(fid, '\\\\ \n')
fprintf(fid,'$e_{aQ}$')
fprintf(fid, '& %3i ', e5aQ(1,:));
fprintf(fid, '\\\\ \n')
fprintf(fid,'$e_{bQ}$')
fprintf(fid, '& %3i ', e5bQ(1,:));
fprintf(fid, '\\\\ \n')

for q = 1:8
  fprintf(fid, ' %3i ', q)
  fprintf(fid, '& %3i ', k(q,:));
  fprintf(fid, '\\\\ \n')
end

fclose(fid)
  






