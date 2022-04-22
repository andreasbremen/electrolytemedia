within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp.Reduced;
function calc_pi "Calculate dimensionless pressure"
  input Modelica.SIunits.Pressure p;

  output Real tau(unit="1") "Dimensionless temperature";
algorithm
  tau := p/IF97.PSTAR1;

end calc_pi;
