within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dtau;
function calc_fhelper "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real f[nLfun-1];

protected
  Real I = calc_I(X);
algorithm

  for i in 1:nLfun-1 loop
  f[i] :=-calc_A_log(T, p)*datafun[i].z^2*sqrt(I);
  end for;

end calc_fhelper;
