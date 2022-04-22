within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes;
function calc_der_g_p
  "Calculates pressure derivative of specific Gibbs free energy of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nLifun] gibbs_dp(unit="J/(kg.Pa)");
protected
  SI.Density[nLifun] di=calc_d_i(T, p);
algorithm
  for i in 1:nLifun loop
    if di[i] > 0 then
      gibbs_dp[i] := 1/di[i];
    end if;
  end for;

end calc_der_g_p;
