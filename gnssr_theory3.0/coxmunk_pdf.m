function gcp_coeff  =  coxmunk_pdf( ws, freq, modelname )
%
%  Updated 5/2016 by JLG to work with SoOp signals - wider range of
%  frequencies.
% 
%  Updated 071309 made switching between models cleaner so less chance of
%  undefined variables if 'modelname' not given.
%  Updated 070109 to include the high-wind-speed model for SJK's 2006 paper
%
%  Reduced Cox and Munk terms in the Gram-Charlier series, based upon the
%  corrections in Wilheit's papers.  WS is measured at the height in C/M
%  (41 meters)
%


if (nargin <= 2)
    modelname = 'coxmunk';
end

if strcmp(modelname,'coxmunk')
  % Cox and munk - the default.
    sigma2c =  0.003 + 1.92e-3 * ws;
    sigma2u =  3.16e-3 * ws;
    c21 = 0.01 - (0.0086 * ws);
    c03 = 0.04 - (0.0330 * ws);
    c40 = 0.40 * ones(size(ws));
    c22 = 0.12 * ones(size(ws));
    c04 = 0.23 * ones(size(ws));
elseif strcmp(modelname,'wu')
    sigma2omni = nan(size(ws));
    lowend = logical(ws <= 7 & ws > 1);
    highend = logical( ws > 7 & ws <=20); 
    sigma2omni(lowend) = 0.01 * ( 3.2 + 1.2 * log( ws (lowend)/7));
    sigma2omni(highend) = 0.01 * ( 3.2 + 6.0 * log( ws (highend)/7));
    sigma2u = sigma2omni/2;
    sigma2c = sigma2omni/2;
    c21 = zeros(size(ws));
    c03 = zeros(size(ws));
    c40 = zeros(size(ws));
    c22 = zeros(size(ws));
    c04 = zeros(size(ws));
elseif strcmp(modelname,'sjk')
    lowend = logical (ws < 3.49);
    highend = logical(ws > 46);
    fu = 6 * log(ws) - 4.0;
    fu(lowend) = ws(lowend);
    fu(highend) = 0.411 * ws(highend);
    sigma2u = ( 0.00316 * fu);
    sigma2c = (0.003 + 0.00192 * fu); 
     c21 = zeros(size(ws));
    c03 = zeros(size(ws));
    c40 = zeros(size(ws));
    c22 = zeros(size(ws));
    c04 = zeros(size(ws));
else
    fprintf(' Invalid model specified \n')
    return  
end


if nargin < 2
  freq = 1.575;  % GHz
end

if (freq < 35)
  if  strcmp(modelname,'sjk')
          sigma2u = 0.45 * sigma2u;
          sigma2c = 0.45 * sigma2c;
  else
     sigma2u = ( 0.3 + 0.02*freq) * sigma2u;
     sigma2c = ( 0.3 + 0.02*freq) * sigma2c;
  end
end



gcp_coeff = [sigma2u, sigma2c, c21, c03, c40, c22, c04];

return

