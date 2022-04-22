within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtau2;
function calc_d2loggamma_dtau2 "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real d2loggamma_dtau2[nLfun - 1];

protected
  Real[nLfun-1] h = calc_hhelper(T,p,X);
  Real[nLfun-1] dh_dtau = calc_dhhelper_dtau(T,p,X);
  Real[nLfun-1] i = calc_ihelper(T,p,X);
  Real[nLfun-1] di_dtau = calc_dihelper_dtau(T,p,X);

algorithm

  for n in 1:nLfun-1 loop
    d2loggamma_dtau2[n] := (dh_dtau[n]*i[n]-di_dtau[n]*h[n])/i[n]^2;
  end for;

end calc_d2loggamma_dtau2;
