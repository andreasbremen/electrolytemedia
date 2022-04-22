within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_cv
  "Calculates specific heat capacity at const volume of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificHeatCapacityAtConstantVolume cv;

protected
  SI.MoleFraction[nGfun] y=calc_Y(X);
algorithm

  cv := Molar.calc_cv(T, d, y)/calc_MM(X);

end calc_cv;
