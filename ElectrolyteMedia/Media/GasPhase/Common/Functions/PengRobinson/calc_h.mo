within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_h
  "Calculates molar enthalpy of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificEnthalpy h;

protected
  SI.MoleFraction[nGfun] y=calc_Y(X);
algorithm

  h := Molar.calc_h(T, d, y)/calc_MM(X);

end calc_h;
