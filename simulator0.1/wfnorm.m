function [wfnorm, nf, delaymean, totalpwr]  = wfnorm( delay, wf )
%
%  Normalizes an ensemble of waveforms (simulated or real) to have
%  a constant average area of unity (in units of delay)
%
%  wf(sample, delay),  delay(sample,delay) or delay(1,delay)
%
nf = mean(wf(:,1));
delaymean = mean(delay,1);
dtau = delaymean(2)-delaymean(1); 
totalpwr = mean(sum((wf-nf),2))*dtau;

wfnorm = (wf - nf)/totalpwr;

end

