within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_d "Unsymmetrical mixing, Pitzer p.124"
  input Real x;
  output Real d[23];
protected
  Real z;
  Real[23] b;
  Integer i;
algorithm

  b := Solute.calc_b(x);
  z := Solute.calc_z(x);

  for j in 1:21 loop
    i :=22 - j;
    d[i] :=b[i + 1] + z*d[i + 1] - d[i + 2];
  end for;

end calc_d;
