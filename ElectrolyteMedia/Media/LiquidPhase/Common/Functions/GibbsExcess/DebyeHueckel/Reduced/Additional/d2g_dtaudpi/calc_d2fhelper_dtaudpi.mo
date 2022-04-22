within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi;
function calc_d2fhelper_dtaudpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real d2f_dtaudpi[nLfun - 1];

protected
  Real I = calc_I(X);
algorithm

  for i in 1:nLfun - 1 loop
    d2f_dtaudpi[i] :=-datafun[i].z^2*sqrt(I)*calc_d2Alog_dtaudpi(T, p);
  end for;

end calc_d2fhelper_dtaudpi;
