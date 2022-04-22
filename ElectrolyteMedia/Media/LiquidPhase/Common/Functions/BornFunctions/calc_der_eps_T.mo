within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_der_eps_T
  "Calculates derivative of dielectric constant w.r.t. T from Johnson1992 from T and rho of water"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real depsdT;
protected
  Real dT_hat = 1 / Born.T_eps;
  Modelica.SIunits.Density rho=IF97_R1_Tp.calc_rho(T, p);
  Real rhoT=IF97_R1_Tp.calc_der_rho_T(T, p);
  Real rho_hat = rho / Born.rho_eps;
  Real T_hat = T / Born.T_eps;
  Real c[:] = {1, Born.a_i[1] / T_hat, Born.a_i[2] / T_hat + Born.a_i[3] + Born.a_i[4] * T_hat, Born.a_i[5] / T_hat + Born.a_i[6] * T_hat + Born.a_i[7] * T_hat ^ 2, Born.a_i[8] / T_hat ^ 2 + Born.a_i[9] / T_hat + Born.a_i[10]};
  Real dcdT[:] = {0, -Born.a_i[1] / T_hat ^ 2 * dT_hat, ((-Born.a_i[2] / T_hat ^ 2) + Born.a_i[4]) * dT_hat, ((-Born.a_i[5] / T_hat ^ 2) + Born.a_i[6] + 2 * Born.a_i[7] * T_hat) * dT_hat, ((-2 * Born.a_i[8] / T_hat ^ 3) - Born.a_i[9] / T_hat ^ 2) * dT_hat};
algorithm
   depsdT := sum(rho_hat ^ (i - 1) * (dcdT[i] + c[i]*(i - 1) * rhoT/rho) for i in 1:5);
annotation(smoothOrder=5);
end calc_der_eps_T;
