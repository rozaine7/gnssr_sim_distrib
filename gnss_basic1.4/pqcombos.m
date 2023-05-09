p = [ 0 1 2 3   0 1 2 3   0 1 2 3   0 1 2 3  ];
q = [ 0 0 0 0   1 1 1 1   2 2 2 2   3 3 3 3  ];
EA = sqrt(2) * exp( (pi/4 + pi/2*p)*i);
EB = sqrt(2) * exp( (pi/4 + pi/2*q)*i);

eai = real(EA);
eaq = imag(EA);
ebi = real(EB);
ebq = imag(EB);

ebar = eai .* eaq .* ebi .* ebq;

