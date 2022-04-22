within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi;
function calc_djhelper_dpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real dj_dpi[nLfun - 1];

protected
   Real[nLfun-1] f = dg_dtau.calc_fhelper(T,p,X);
   Real[nLfun-1] g = dg_dtau.calc_ghelper(T,p,X);
   Real[nLfun-1] df_dtau = dg_dtau.calc_dfhelper_dtau(T,p,X);
   Real[nLfun-1] df_dpi = dg_dpi.calc_dfhelper_dpi(T,p,X);
   Real[nLfun-1] dghelper_dtau = dg_dtau.calc_dghelper_dtau(T,p,X);
   Real[nLfun-1] dghelper_dpi = dg_dpi.calc_dghelper_dpi(T,p,X);
   Real[nLfun-1] d2f_dtaudpi = calc_d2fhelper_dtaudpi(T,p,X);
   Real[nLfun-1] d2g_dtaudpi_ = calc_d2ghelper_dtaudpi(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    dj_dpi[i] :=d2f_dtaudpi[i]*g[i]+df_dtau[i]*dghelper_dpi[i]-(df_dpi[i]*dghelper_dtau[i]+f[i]*d2g_dtaudpi_[i]);
  end for;

end calc_djhelper_dpi;
