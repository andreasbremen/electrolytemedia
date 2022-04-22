within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_2der_eps_pT
  "Calculates second derivative of dielectric constant w.r.t. T and p from Johnson1992 from T and rho of water"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real epspT;
protected
  parameter Real dT_hat = 1 / Born.T_eps;
  Modelica.SIunits.Density rho = IF97_R1_Tp.calc_rho(T,p);
  Real rhop=IF97_R1_Tp.calc_der_rho_p(T, p);
  Real rhoT=IF97_R1_Tp.calc_der_rho_T(T, p);
  Real rhopT=IF97_R1_Tp.calc_2der_rho_pT(T, p);
  Real rho_hat = rho / Born.rho_eps;
  Real T_hat = T / Born.T_eps;
  Real c[:] = {1, Born.a_i[1] / T_hat, Born.a_i[2] / T_hat + Born.a_i[3] + Born.a_i[4] * T_hat, Born.a_i[5] / T_hat + Born.a_i[6] * T_hat + Born.a_i[7] * T_hat ^ 2, Born.a_i[8] / T_hat ^ 2 + Born.a_i[9] / T_hat + Born.a_i[10]};
  Real dcdT[:] = {0, -Born.a_i[1] / T_hat ^ 2 * dT_hat, ((-Born.a_i[2] / T_hat ^ 2) + Born.a_i[4]) * dT_hat, ((-Born.a_i[5] / T_hat ^ 2) + Born.a_i[6] + 2 * Born.a_i[7] * T_hat) * dT_hat, ((-2 * Born.a_i[8] / T_hat ^ 3) - Born.a_i[9] / T_hat ^ 2) * dT_hat};
algorithm
   epspT := sum((i-1)*rho_hat ^ (i - 1) * (c[i]/rho*rhopT - c[i]/rho^2*rhoT*rhop + rhop/rho*(dcdT[i]+c[i]*(i-1)*rhoT/rho)) for i in 1:5);
annotation(smoothOrder=5);
end calc_2der_eps_pT;
