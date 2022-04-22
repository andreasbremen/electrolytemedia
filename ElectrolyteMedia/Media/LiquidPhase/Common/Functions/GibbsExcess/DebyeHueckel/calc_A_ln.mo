within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel;
function calc_A_ln
  "Calculates Debye Hückel parameter A_Phi in ln basis from Pitzer1991"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real A;
protected
  SI.Density rho = IF97_R1_Tp.calc_rho(T,p);
  Real eps=BornFunctions.calc_eps(T,p);
  Real e=Modelica.Constants.F/Modelica.Constants.N_A
    "electronic charge 1.602176621e-19";
algorithm
  A := 1/3*(2*Modelica.Constants.pi*Modelica.Constants.N_A*rho)^0.5*(e^2/(4*
    Modelica.Constants.pi*eps*Modelica.Constants.epsilon_0*Modelica.Constants.k*
    T))^1.5;
end calc_A_ln;
