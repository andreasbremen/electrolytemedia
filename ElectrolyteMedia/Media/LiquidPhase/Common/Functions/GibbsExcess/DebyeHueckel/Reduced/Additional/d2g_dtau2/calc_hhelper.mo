within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtau2;
function calc_hhelper "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real h[nLfun - 1];

protected
  Real[nLfun - 1] dfhelper_dtau = dg_dtau.calc_dfhelper_dtau(T,p,X);
  Real[nLfun - 1] f = dg_dtau.calc_fhelper(T,p,X);
  Real[nLfun - 1] dghelper_dtau = dg_dtau.calc_dghelper_dtau(T,p,X);
  Real[nLfun - 1] g = dg_dtau.calc_ghelper(T,p,X);

algorithm

   for i in 1:nLfun-1 loop
     h[i] := dfhelper_dtau[i]*g[i]-dghelper_dtau[i]*f[i];
   end for;

end calc_hhelper;
