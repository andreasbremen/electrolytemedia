within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_J_dx "Unsymmetrical mixing, Pitzer p.124"
  input Real x;
  output Real J_dx;
protected
  Real[23] d;
  Real z_dx;
algorithm
  z_dx := Solute.calc_z_dx(x);
  d := Solute.calc_d(x);

  J_dx :=0.25 + 0.5* z_dx *(d[1] - d[3]);

end calc_J_dx;
