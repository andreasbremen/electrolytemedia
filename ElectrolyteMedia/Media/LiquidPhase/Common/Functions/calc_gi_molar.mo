within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_gi_molar
  "Calculates partial Gibbs free energies of solutes and solvent at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.MolarEnergy[nLfun] gi;

algorithm
  gi[1:nLfun-1] :=Solutes.calc_g_i(T, p) .* datafun[:].MM;
  gi[nLfun] :=IF97_R1_Tp.calc_g(T, p)*IF97.MH2O;

  annotation(smoothOrder=5);
end calc_gi_molar;
