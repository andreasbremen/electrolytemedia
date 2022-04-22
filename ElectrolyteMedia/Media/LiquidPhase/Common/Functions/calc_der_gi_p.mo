within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_der_gi_p
  "Calculates pressure derivatives of partial Gibbs free energies of solutes and solvent at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificEnergy[nLifun+1] gi;

algorithm
  gi[1:nLifun] :=Solutes.calc_der_g_p(T, p);
  gi[nLifun+1] :=IF97_R1_Tp.calc_der_g_p(T, p);

end calc_der_gi_p;
