within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced;
function calc_2der_g_taupi
  "Calculates 2nd derivative w.r.t. tau,pi of excess reduced Gibbs free energy at infinite dilution of aqueous species"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real d2g_dtaudpi[nLfun];

protected
  Real ai[nLfun - 1]=calc_activity_i(
              T,
              p,
              X);
  Real dai_dpi[nLfun - 1]=Additional.dg_dpi.calc_dai_dpi(
      T,
      p,
      X);
  Real dai_dtau[nLfun - 1]=Additional.dg_dtau.calc_dai_dtau(
      T,
      p,
      X);
  Real d2ai_dtaudpi[nLfun - 1]=Additional.d2g_dtaudpi.calc_d2ai_dtaudpi(
      T,
      p,
      X);

algorithm

  for i in 1:nLfun-1 loop
    if ai[i] > 0 then
      d2g_dtaudpi[i] := -ai[i]^(-2)*dai_dpi[i]*dai_dtau[i]+ai[i]^(-1)*d2ai_dtaudpi[i];
    end if;
  end for;

  d2g_dtaudpi[nLfun] :=0;

end calc_2der_g_taupi;
