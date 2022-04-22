within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_theta_0 "Calculates Einstein temperature at T_0"

  output Real theta_0[nSfun];
algorithm
  theta_0 := 10636 ./ (datafun[:].S_ref ./ datafun[:].n + 6.44 * ones(nSfun));
end calc_theta_0;
