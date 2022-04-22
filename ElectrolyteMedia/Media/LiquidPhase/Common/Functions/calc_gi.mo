within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_gi
  "Calculates partial Gibbs free energies of solutes and solvent at T and p at infinite dilution"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificEnergy[nLfun] gi;

algorithm
  gi[1:nLfun-1] :=Solutes.calc_g_i(T, p);
  gi[nLfun] :=IF97_R1_Tp.calc_g(T, p);

  annotation(smoothOrder=5);
end calc_gi;
