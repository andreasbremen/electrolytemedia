within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel;
function calc_der_g_p
  "Calculates pressure derivative of specific Gibbs free energy of solutes at inifinite dilution based on HKF model"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Real[nLfun] gibbs_dp(unit="J/(kg.Pa)");
protected
  SI.Density[nLfun] di = calc_d(T,p,X);
algorithm
  for i in 1:nLfun loop
    if di[i] > 0 then
      gibbs_dp[i] := 1/di[i];
    end if;
  end for;

end calc_der_g_p;
