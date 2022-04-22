within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel;
function calc_A_log
  "Calculates Debye Hückel parameters A_gamma in (dm3/mol)^0.5 and B_gamma in (m/mol)^0.5"

  input SI.Temperature T;
  input SI.Pressure p;
  output Real A_gamma_log;
protected
  Real A_Phi_ln = calc_A_ln(T,p);
algorithm
  A_gamma_log := A_Phi_ln * 3 / log(10);
end calc_A_log;
