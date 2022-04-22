within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp.Reduced;
function calc_3der_g_taupipi
  "Third derivative of Gibbs free energy w.r.t. tau,pi and pi"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real gpipipi;
protected
  Real pi=Reduced.calc_pi(p);
algorithm
  gpipipi :=p/pi*calc_3der_g_taupip(T, p);
end calc_3der_g_taupipi;
