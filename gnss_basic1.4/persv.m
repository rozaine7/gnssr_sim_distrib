function [rt, pr, cyc, phase, slp_dtct, cn0] = persv( gps_mat )
%
% function to extract the measurements per satellite PRN from the data
% given in Misra and Enge's book.  One column per SV - the flag of -999
% indicates that no data was received for that SV.
%
%nrows = min(size(gps_mat,1), nreq);
nrows = size(gps_mat,1);
rt(1,1) = gps_mat(1,1);
for k=2:nrows
  if (gps_mat(k,1) > rt(end))
    rt = [rt; gps_mat(k,1)];
  end
end

nobs = size(rt,1);
pr = -999 * ones(nobs, 32);
cyc = -999 * ones(nobs, 32);
phase = -999 * ones(nobs, 32);
slp_dtct = -999 * ones(nobs, 32);
cn0 = -999 * ones(nobs, 32);

for k=1:nobs
  svset = gps_mat(logical(gps_mat(:,1) == rt(k)),2)';
  prset = gps_mat(logical(gps_mat(:,1) == rt(k)),3)';
  pr(k, svset) = prset;
  cycset = gps_mat(logical(gps_mat(:,1) == rt(k)),4)';
  cyc(k, svset) = cycset;
  phset = gps_mat(logical(gps_mat(:,1) == rt(k)),5)';
  phase(k, svset) = phset;
  sdset = gps_mat(logical(gps_mat(:,1) == rt(k)),6)';
  slp_dtct(k, svset) = sdset;
  cn0set = gps_mat(logical(gps_mat(:,1) == rt(k)),7)';
  cn0(k, svset) = cn0set;
end

return