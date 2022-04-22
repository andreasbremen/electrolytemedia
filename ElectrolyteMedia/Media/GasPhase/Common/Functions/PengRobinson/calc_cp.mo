within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_cp
  "Calculates specificheat capacity of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificHeatCapacityAtConstantPressure cp;

protected
  SI.MoleFraction[nGfun] y=calc_Y(X);
algorithm

  cp := Molar.calc_cp(T, d, y)/calc_MM(X);

end calc_cp;
