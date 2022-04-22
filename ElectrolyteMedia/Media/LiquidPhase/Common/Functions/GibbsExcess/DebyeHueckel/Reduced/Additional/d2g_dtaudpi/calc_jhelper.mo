within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi;
function calc_jhelper "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real j[nLfun - 1];

protected
  Real[nLfun-1] f = dg_dtau.calc_fhelper(T,p,X);
  Real[nLfun-1] g = dg_dtau.calc_ghelper(T,p,X);
  Real[nLfun-1] f_= dg_dtau.calc_dfhelper_dtau(T,p,X);
  Real[nLfun-1] g_ = dg_dtau.calc_dghelper_dtau(T,p,X);
algorithm

  for i in 1:nLfun-1 loop
    j[i] :=f_[i]*g[i] - f[i]*g_[i];
  end for;

end calc_jhelper;
