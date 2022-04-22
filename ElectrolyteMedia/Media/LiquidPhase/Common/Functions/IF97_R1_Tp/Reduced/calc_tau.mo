within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp.Reduced;
function calc_tau "Calculate dimensionless temperature"
  input Modelica.SIunits.Temperature T;

  output Real tau(unit="1") "Dimensionless temperature";
algorithm
  tau := IF97.TSTAR1/T;

end calc_tau;
