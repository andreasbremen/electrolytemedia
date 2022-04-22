within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dpi2;
function calc_hhelper "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real h[nLfun - 1];

protected
  Real[nLfun - 1] dfhelper_dpi = dg_dpi.calc_dfhelper_dpi(T,p,X);
  Real[nLfun - 1] f = dg_dpi.calc_fhelper(T,p,X);
  Real[nLfun - 1] dghelper_dpi = dg_dpi.calc_dghelper_dpi(T,p,X);
  Real[nLfun - 1] g = dg_dpi.calc_ghelper(T,p,X);

algorithm

   for i in 1:nLfun-1 loop
     h[i] := dfhelper_dpi[i]*g[i]-dghelper_dpi[i]*f[i];
   end for;

end calc_hhelper;
