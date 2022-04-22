within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dpi;
function calc_dloggamma_dpi "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real dloggamma_dpi[nLfun-1];

protected
  Real f[nLfun - 1]=calc_fhelper(T,p,X);
  Real df_dpi[nLfun - 1]=calc_dfhelper_dpi(T,p,X);
  Real g[nLfun - 1]=calc_ghelper(T,p,X);
  Real dg_dpi[nLfun - 1]=calc_dghelper_dpi(T,p,X);

algorithm

  for i in 1:nLfun - 1 loop
    dloggamma_dpi[i] := (df_dpi[i]*g[i] -
    dg_dpi[i]*f[i])/g[i]^2;
  end for;

end calc_dloggamma_dpi;
