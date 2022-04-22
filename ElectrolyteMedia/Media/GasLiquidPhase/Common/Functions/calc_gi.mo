within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_gi
  "calculates partial Gibbs free energies of solutes and solvent at T and p"
  input SI.Temperature T;
  input SI.Pressure p;

  output SI.SpecificEnergy[nGfunfun+nLfunfun] gi;
algorithm
   gi[1:nGfunfun] :=GasFunctions.calc_g_i(T);
   gi[1+nGfunfun:nGfunfun+nLfunfun] :=LiquidFunctions.calc_gi(T,p);

  annotation(Inline=true,smoothOrder=5);
end calc_gi;
