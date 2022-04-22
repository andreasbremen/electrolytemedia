within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dpi2;
function calc_d2loggamma_dpi2 "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real d2loggamma_dpi2[nLfun - 1];

protected
  Real[nLfun-1] h = calc_hhelper(T,p,X);
  Real[nLfun-1] dh_dpi = calc_dhhelper_dpi(T,p,X);
  Real[nLfun-1] i = calc_ihelper(T,p,X);
  Real[nLfun-1] di_dpi = calc_dihelper_dpi(T,p,X);

algorithm

  for n in 1:nLfun-1 loop
    d2loggamma_dpi2[n] := (dh_dpi[n]*i[n]-di_dpi[n]*h[n])/i[n]^2;
  end for;
end calc_d2loggamma_dpi2;
