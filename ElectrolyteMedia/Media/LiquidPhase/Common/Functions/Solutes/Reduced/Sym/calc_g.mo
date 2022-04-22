within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced.Sym;
function calc_g
  "Calculates reduced Gibbs free energy for asymmetrical mixing"

  input SI.Temperature T;
  input SI.MassFraction[nLfun] X;
  output Real g;

protected
  Real [nLfun] Y = calc_Y(X);

algorithm

  g:= (1 - Y[nLfun])*log(1/IF97.MH2O);

end calc_g;
