within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_der_eps_TT
  "Calculates second derivative of dielectric constant w.r.t. T from Johnson1992 from T and rho of water"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real d2epsdT;
protected
  Real dT_hat = 1 / Born.T_eps;
  Modelica.SIunits.Density rho=IF97_R1_Tp.calc_rho(T, p);
  Real rhoT=IF97_R1_Tp.calc_der_rho_T(T, p);
  Real rhoTT=IF97_R1_Tp.calc_2der_rho_TT(T, p);
  Real rho_hat = rho / Born.rho_eps;
  Real T_hat = T / Born.T_eps;
  Real c[:] = {1, Born.a_i[1] / T_hat, Born.a_i[2] / T_hat + Born.a_i[3] + Born.a_i[4] * T_hat, Born.a_i[5] / T_hat + Born.a_i[6] * T_hat + Born.a_i[7] * T_hat ^ 2, Born.a_i[8] / T_hat ^ 2 + Born.a_i[9] / T_hat + Born.a_i[10]};
  Real dcdT[:] = {0, -Born.a_i[1] / T_hat ^ 2 * dT_hat, ((-Born.a_i[2] / T_hat ^ 2) + Born.a_i[4]) * dT_hat, ((-Born.a_i[5] / T_hat ^ 2) + Born.a_i[6] + 2 * Born.a_i[7] * T_hat) * dT_hat, ((-2 * Born.a_i[8] / T_hat ^ 3) - Born.a_i[9] / T_hat ^ 2) * dT_hat};
  Real d2cdT[:] = {0, 2 * Born.a_i[1] / T_hat ^ 3 * dT_hat ^ 2, 2 * Born.a_i[2] / T_hat ^ 3 * dT_hat ^ 2, 2 * Born.a_i[5] / T_hat ^ 3 * dT_hat ^ 2 + 2 * Born.a_i[7] * dT_hat ^ 2, (6 * Born.a_i[8] / T_hat ^ 4 + 2 * Born.a_i[9] / T_hat ^ 3) * dT_hat ^ 2};
algorithm
  d2epsdT := sum(rho_hat ^ (i - 1) * (d2cdT[i] + (i - 1)/rho * (2*dcdT[i]*rhoT + c[i]*rhoTT + c[i]*(i-2)*rhoT^2/rho)) for i in 1:5);
annotation(smoothOrder=5);
end calc_der_eps_TT;
