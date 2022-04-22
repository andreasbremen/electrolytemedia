within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_J "Unsymmetrical mixing, Pitzer p.124"
  input Real x;
  output Real J;
protected
  Real[23] b;
algorithm
  b := Solute.calc_b(x);

  J :=0.25*x - 1 + 0.5*(b[1] - b[3]);

end calc_J;
