within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_MM "Returns molar mass of mixture"

  input SI.MassFraction X[nGfun];
  output SI.MolarMass MM;

protected
  SI.MoleFraction[nGfun] y=calc_Y(X);

algorithm
  MM := Molar.calc_MM(y);
end calc_MM;
