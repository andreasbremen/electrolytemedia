within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dtau;
function calc_dfhelper_dtau "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real df_dT[nLfun - 1];

protected
  Real I = calc_I(X);
algorithm

  for i in 1:nLfun - 1 loop
    df_dT[i] :=-datafun[i].z^2*sqrt(I)*calc_dAlog_dtau(T, p);
  end for;

end calc_dfhelper_dtau;
