within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel;
function calc_der_g_T
  "Calculates temperature derivative of specific Gibbs free energy of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Real[nLfun] gibbs_dT(unit="J/(kg.K)");
algorithm
  gibbs_dT := - calc_s(T,p,X);

end calc_der_g_T;
