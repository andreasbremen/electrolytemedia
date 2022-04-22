within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_der_eps_pp
  "Calculates derivative of dielectric constant w.r.t. p from Johnson1992 from T and rho of water"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real epspp;
protected
  Modelica.SIunits.Density rho = IF97_R1_Tp.calc_rho(T,p);
  Real rhop=IF97_R1_Tp.calc_der_rho_p(T, p);
  Real rhopp=IF97_R1_Tp.calc_2der_rho_pp(T, p);
  Real rho_hat = rho / Born.rho_eps;
  Real T_hat = T / Born.T_eps;
  Real c[:] = {1, Born.a_i[1] / T_hat, Born.a_i[2] / T_hat + Born.a_i[3] + Born.a_i[4] * T_hat, Born.a_i[5] / T_hat + Born.a_i[6] * T_hat + Born.a_i[7] * T_hat ^ 2, Born.a_i[8] / T_hat ^ 2 + Born.a_i[9] / T_hat + Born.a_i[10]};
algorithm
  epspp := sum((i - 1) * c[i] * rho_hat ^ (i - 1)*(1/rho*rhopp + (i-2)/rho^2*rhop^2) for i in 1:5);
annotation(smoothOrder=5);
end calc_der_eps_pp;
