within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes;
function calc_g_i
  "Calculates specific Gibbs free energy of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  output SI.SpecificGibbsFreeEnergy[nLifun] gibbs;
algorithm
  gibbs :=Reduced.calc_g(T,p).*datafun.R*T;

  annotation(smoothOrder=5);
end calc_g_i;
