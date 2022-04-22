within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_eps
  "Calculates dielectric constant from Johnson1992 from T and rho of water"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real eps;
protected
  SI.Density rho=IF97_R1_Tp.calc_rho(T, p);
  Real rho_hat = rho / Born.rho_eps;
  Real T_hat = T / Born.T_eps;
  Real c[:] = {1, Born.a_i[1] / T_hat, Born.a_i[2] / T_hat + Born.a_i[3] + Born.a_i[4] * T_hat, Born.a_i[5] / T_hat + Born.a_i[6] * T_hat + Born.a_i[7] * T_hat ^ 2, Born.a_i[8] / T_hat ^ 2 + Born.a_i[9] / T_hat + Born.a_i[10]};
algorithm
  eps := sum(c[i] * rho_hat ^ (i - 1) for i in 1:5);
annotation(smoothOrder=20);
end calc_eps;
