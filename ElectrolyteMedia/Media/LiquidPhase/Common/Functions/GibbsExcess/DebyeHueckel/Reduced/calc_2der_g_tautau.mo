within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced;
function calc_2der_g_tautau
  "Calculates 2nd derivative w.r.t. tau of excess reduced Gibbs free energy at infinite dilution of aqueous species"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real d2g_dtau2[nLfun];

protected
  Real[nLfun - 1] ai=calc_activity_i(
              T,
              p,
              X);
  Real[nLfun - 1] dai_dtau=Additional.dg_dtau.calc_dai_dtau(
      T,
      p,
      X);
  Real[nLfun - 1] d2ai_dtau2=Additional.d2g_dtau2.calc_d2ai_dtau2(
      T,
      p,
      X);

algorithm

  for i in 1:nLfun-1 loop
    if ai[i] > 0 then
      d2g_dtau2[i] := -ai[i]^(-2)*dai_dtau[i]^2+ai[i]^(-1)*d2ai_dtau2[i];
    end if;
  end for;
  d2g_dtau2[nLfun] := 0;

end calc_2der_g_tautau;
