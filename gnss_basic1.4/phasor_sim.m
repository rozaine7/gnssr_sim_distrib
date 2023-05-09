function [Ihist, Qhist, thist] = phasor_sim( data_set, sample_rate, ...
			fcarrier, block_size, pnmatrix, pnflag);
nsamples = size(data_set, 1);
nblocks = floor(nsamples / block_size);
Ihist = zeros(nblocks,1);
Qhist = zeros(nblocks,1);

[timetags, baseband_local, cid] = bpsk( sample_rate, nsamples, pnmatrix, ...
					pnflag); 

[Ilocal, Qlocal] = carrier( sample_rate, baseband_local, fcarrier, 0);

for i=1:nblocks
  Ihist(i,1)  = td_corr( data_set, Ilocal, (i-1)*block_size+1, 0, ...
			 block_size);
  Qhist(i,1)  = td_corr( data_set, Qlocal, (i-1)*block_size+1, 0, ...
			   block_size); 
  thist(i,1) = timetags(i*block_size);
end