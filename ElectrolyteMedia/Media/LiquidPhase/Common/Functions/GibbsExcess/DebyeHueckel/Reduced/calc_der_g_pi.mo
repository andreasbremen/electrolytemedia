within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced;
function calc_der_g_pi
  "Calculates derivative w.r.t. pi of excess reduced Gibbs free energy at infinite dilution of aqueous species"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real dg_dpi[nLfun];

protected
  Real ai[nLfun - 1]=calc_activity_i(
              T,
              p,
              X);
  Real dai_dpi[nLfun - 1]=Additional.dg_dpi.calc_dai_dpi(
      T,
      p,
      X);

algorithm
  for i in 1:nLfun-1 loop
    if ai[i] > 0 then
      dg_dpi[i] := 1/(ai[i]) * dai_dpi[i];
    end if;
  end for;
  dg_dpi[nLfun] := 0;

end calc_der_g_pi;
