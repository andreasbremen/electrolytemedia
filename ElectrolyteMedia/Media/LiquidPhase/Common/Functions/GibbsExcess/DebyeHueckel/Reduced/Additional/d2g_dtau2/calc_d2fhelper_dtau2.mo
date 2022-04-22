within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtau2;
function calc_d2fhelper_dtau2 "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real d2f_dT2[nLfun - 1];

protected
  Real I = calc_I(X);
algorithm

  for i in 1:nLfun - 1 loop
    d2f_dT2[i] :=-datafun[i].z^2*sqrt(I)*calc_d2Alog_dtau2(T, p);
  end for;

end calc_d2fhelper_dtau2;
