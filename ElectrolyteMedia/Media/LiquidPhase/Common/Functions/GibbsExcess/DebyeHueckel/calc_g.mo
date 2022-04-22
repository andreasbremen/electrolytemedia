within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel;
function calc_g "Calculates specific Gibbs free energy with Debye Hückel model from temperature, pressure and mass fraction"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output SI.SpecificEnergy[nLfun] g;

protected
  Real[nLfun] gred = Reduced.calc_g(T,p,X);
  SI.SpecificEntropy[nLfun] R = calc_RX();

algorithm

  g :=T*R .* gred;

end calc_g;
