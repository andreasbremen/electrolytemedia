within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_der_gi_T
  "Calculates temperature derivatives of partial Gibbs free energies of solutes and solvent at T and p at infinite dilution"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificEnergy[nLifun+1] gi;

algorithm
  gi[1:nLifun] :=Solutes.calc_der_g_T(T, p);
  gi[nLifun+1] :=IF97_R1_Tp.calc_der_g_T(T, p);

end calc_der_gi_T;
