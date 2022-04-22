within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dtau;
function calc_ghelper "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real g[nLfun-1];

protected
  Real I = calc_I(X);
algorithm

  for i in 1:nLfun-1 loop
  g[i] :=1 + calc_B_log(T, p)*1e-10*3.72*sqrt(I);//datafun[i].a_0
  end for;

end calc_ghelper;
