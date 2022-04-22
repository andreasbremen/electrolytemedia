within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dpi2;
function calc_ihelper "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real i[nLfun - 1];

protected
  Real g[nLfun-1] = dg_dpi.calc_ghelper(T,p,X);

algorithm

  for n in 1:nLfun-1 loop
    i[n] := g[n]^2;
  end for;

end calc_ihelper;
