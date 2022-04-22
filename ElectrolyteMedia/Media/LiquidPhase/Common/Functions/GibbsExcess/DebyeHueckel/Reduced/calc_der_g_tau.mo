within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced;
function calc_der_g_tau
  "Calculates derivative w.r.t. tau of excess reduced Gibbs free energy at infinite dilution of aqueous species"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real dmu_dT[nLfun];

protected
  Real ai[nLfun - 1]=calc_activity_i(
              T,
              p,
              X);
  Real dai_dT[nLfun - 1]=Additional.dg_dtau.calc_dai_dtau(
      T,
      p,
      X);

algorithm
  for i in 1:nLfun-1 loop
    if ai[i] > 0 then
      dmu_dT[i] := 1/(ai[i]) * dai_dT[i];
    end if;
  end for;
  dmu_dT[nLfun] := 0;

end calc_der_g_tau;
