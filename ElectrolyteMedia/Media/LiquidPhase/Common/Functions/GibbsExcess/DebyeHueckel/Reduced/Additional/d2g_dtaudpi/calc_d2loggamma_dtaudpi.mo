within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi;
function calc_d2loggamma_dtaudpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real d2loggamma_dtaudpi[nLfun - 1];

protected
  Real[nLfun-1] j = calc_jhelper(T,p,X);
  Real[nLfun-1] dj_dpi = calc_djhelper_dpi(T,p,X);
  Real[nLfun-1] k = calc_khelper(T,p,X);
  Real[nLfun-1] dk_dpi = calc_dkhelper_dpi(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    d2loggamma_dtaudpi[i] :=(dj_dpi[i]*k[i]-j[i]*dk_dpi[i])/k[i]^2;
  end for;

end calc_d2loggamma_dtaudpi;
