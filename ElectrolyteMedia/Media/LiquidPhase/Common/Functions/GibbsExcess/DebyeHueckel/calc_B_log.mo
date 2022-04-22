within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel;
function calc_B_log
  "Calculates Debye Hückel parameters B_gamma in (m/mol)^0.5"

  extends Modelica.Icons.Function;
  input SI.Temperature T;
  input SI.Pressure p;
  output Real B;
protected
  SI.Density rho = IF97_R1_Tp.calc_rho(T,p);
  Real eps = BornFunctions.calc_eps(T,p);
  Real e=Modelica.Constants.F/Modelica.Constants.N_A
    "electronic charge 1.602176621e-19";
algorithm
  B := (2*e^2*rho*Modelica.Constants.N_A/(Modelica.Constants.epsilon_0*eps*
    Modelica.Constants.k*T))^0.5;
end calc_B_log;
