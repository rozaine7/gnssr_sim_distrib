function [pnmatrix, G1, G2] = pngen(prn);
%
% function to generate the PN codes in binary (0, 1)
%
% prn is a matrix of satellite numbers (PRN's actually)
% pnmatrix is the output array of dimension 1023 X size(prn)
%
% G1 is the maximal length shift prn code from shift register G1.
% G2 is the full state of the 10 stage shift register G2 (1023 X 10)
%
% Generate the G1 codes.
%
G1 = zeros(1023,1);
G2 = zeros(1023,10);
G1_register = [1 1 1 1 1 1 1 1 1 1];
G2_register = [1 1 1 1 1 1 1 1 1 1];
for G1_counter = 1:1023
  G1(G1_counter,1) = G1_register(10);
  feedback1 = xor(G1_register(3), G1_register(10));
  G1_register = [feedback1, G1_register(1:9)];
end
%
% Generate the G2 codes
%
for G2_counter = 1:1023
  G2(G2_counter,:) = G2_register;
  feedback2 = mod(sum(G2_register([2,3,6,8,9,10])),2);
  G2_register = [feedback2, G2_register(1:9)];
end
%
% Generate each of the 32 PRNS
%

taps1 = [2, 3, 4, 5, 1, 2, 1, 2, 3, 2, 3, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, ...
	 6, 1, 4, 5, 6, 7, 8, 1, 2, 3, 4];
taps2 = [6, 7, 8, 9, 9, 10, 8, 9, 10, 3, 4, 6, 7, 8, 9, 10, 4, 5, 6, 7, ...
	 8, 9, 3, 6, 7, 8, 9, 10, 6, 7, 8, 9];

pnmatrix = mod( (G1*ones(1,32) + G2(:,taps1) + G2(:,taps2)),2);

return
