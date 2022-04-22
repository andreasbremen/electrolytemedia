within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Debye;
function calc_g
  "Calculates mixture properties of Gibbs free energy with Debye limiting law"
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
