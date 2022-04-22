within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_theta_T "Calculates dimensionless temperature"

  input SI.Temperature T;

  output Real u[nSfun];
protected
  parameter Real theta_0[nSfun]=Molar.calc_theta_0();
algorithm
  u := theta_0/ T;
  annotation(smoothOrder=2);
end calc_theta_T;
