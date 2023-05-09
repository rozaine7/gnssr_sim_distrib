function r = gauss_corr(xstate, acf, lag)
  var = xstate(1);
  corrtime = xstate(2);
  r = var * exp(-(lag./corrtime).^2) - acf;
end

