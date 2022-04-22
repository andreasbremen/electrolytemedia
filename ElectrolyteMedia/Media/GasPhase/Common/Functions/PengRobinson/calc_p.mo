within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_p
  "Calculates pressure of gas mixture with Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output SI.Pressure p;
protected
  SI.MoleFraction[nGfun] y=calc_Y(X);
algorithm
  p := Molar.calc_p(T,d,y);
  annotation(smoothOrder=5);
end calc_p;
