within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_g_i
  "Calculates specific Gibbs free energy at temperature T with corresponding EOS"
  input SI.Temperature T;

  output SI.SpecificEnergy[nGfun] g;
algorithm
  g :=IdealGas.calc_g_i(T);
annotation(smoothOrder=5);
end calc_g_i;
