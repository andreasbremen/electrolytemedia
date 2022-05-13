within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_gi
  "calculates partial Gibbs free energies of solutes and solvent at T and p"
  input SI.Temperature T;
  input SI.Pressure p;

  output SI.SpecificEnergy[nSfunfun+nLfunfun] gi;
algorithm
   gi[1:nSfunfun] :=SolidFunctions.calc_g_i(T,p);
   gi[1+nSfunfun:nSfunfun+nLfunfun] :=LiquidFunctions.calc_gi(T,p);

  annotation(smoothOrder=5);
end calc_gi;
