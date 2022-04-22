within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_u
  "Calculates internal energy of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificEnergy u;

protected
  SI.MoleFraction[nGfun] y=calc_Y(X);
algorithm

  u := Molar.calc_u(T, d, y)/calc_MM(X);
end calc_u;
