within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_Z_dpdT
  "Calculates second derivative w.r.t. p and T of Born function Z"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real ZpT;
protected
  Real eps=calc_eps(T, p);
  Real epsT=calc_der_eps_T(T, p);
  Real epsp=calc_der_eps_p(T, p);
  Real epspT=calc_2der_eps_pT(T, p);
algorithm
  ZpT :=(eps*epspT - 2*epsT*epsp)/eps^3;
  annotation(smoothOrder=5);
end calc_Z_dpdT;
