within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi;
function calc_khelper "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real k[nLfun - 1];

protected
  Real[nLfun-1] g=dg_dpi.calc_ghelper(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    k[i] :=g[i]^2;
  end for;

end calc_khelper;
