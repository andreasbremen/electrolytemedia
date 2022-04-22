within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes;
function calc_der_g_T
  "Calculates temperature derivative of specific Gibbs free energy of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nLifun] gibbs_dT(unit="J/(kg.K)");
algorithm
  gibbs_dT :=-calc_s_i(T, p);

end calc_der_g_T;
