within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dtau;
function calc_dloggamma_dtau "Helper function to calculate Gibbs derivative"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real dloggamma_dT[nLfun-1];

protected
  Real f[nLfun - 1]=calc_fhelper(T,p,X);
  Real df_dT[nLfun - 1]=calc_dfhelper_dtau(T,p,X);
  Real g[nLfun - 1]=calc_ghelper(T,p,X);
  Real dg_dT[nLfun - 1]=calc_dghelper_dtau(T,p,X);

algorithm

  for i in 1:nLfun - 1 loop
    dloggamma_dT[i] := (df_dT[i]*g[i] -
    dg_dT[i]*f[i])/g[i]^2;
  end for;

end calc_dloggamma_dtau;
