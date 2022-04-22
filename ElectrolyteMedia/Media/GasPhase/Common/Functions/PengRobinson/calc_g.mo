within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_g
  "Calculates specific Gibbs free energy of gas phase with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output SI.SpecificEnergy g;

protected
  SI.MoleFraction[nGfun] y=calc_Y(X);
algorithm

  g := Molar.calc_g(T, d, y)/calc_MM(X);
end calc_g;
