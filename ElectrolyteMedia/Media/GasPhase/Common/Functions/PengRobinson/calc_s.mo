within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_s "Calculates entropy of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificEntropy s;

protected
  SI.MoleFraction[nGfun] y=calc_Y(X);
algorithm

  s := Molar.calc_s(T, d, y)/calc_MM(X);

end calc_s;
