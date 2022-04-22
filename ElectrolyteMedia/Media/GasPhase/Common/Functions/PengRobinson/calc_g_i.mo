within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_g_i "Calculates specific Gibbs free energy at temperature T"

  input SI.Temperature T "Temperature";
  output SI.SpecificGibbsFreeEnergy[nGfun] g "Specific Gibbs free energy at temperature T";
algorithm
  g :=IdealGas.Molar.calc_g_i(T) ./ datafun.MM;
  annotation(smoothOrder=5);
end calc_g_i;
