within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_g "Calculates specific Gibbs free energy of a mixture of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfun] X;
  output SI.SpecificEnergy g_mix;

protected
  SI.SpecificEnergy g[nSfun]=calc_g_i(T, p);

algorithm

  g_mix :=g*X;

end calc_g;
