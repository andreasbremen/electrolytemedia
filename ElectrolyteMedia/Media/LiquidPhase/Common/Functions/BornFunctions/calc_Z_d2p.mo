within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.BornFunctions;
function calc_Z_d2p "Calculates second p derivative of Born function Z"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  output Real Z;
protected
  Real eps=calc_eps(T, p);
  Real epsp=calc_der_eps_p(T, p);
  Real epspp=calc_der_eps_pp(T, p);
algorithm
  Z := (eps*epspp - 2*epsp^2)/eps^3;
  annotation(smoothOrder=5);
end calc_Z_d2p;
