within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_g
  "calculates specific Gibbs free energy of gas, solid mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output SI.SpecificGibbsFreeEnergy g;

protected
  SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
  SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
  SI.SpecificGibbsFreeEnergy gs = SolidFunctions.calc_g(T,p,Xs);
  SI.SpecificGibbsFreeEnergy gl = LiquidFunctions.calc_g(T,p,Xl);
algorithm

  g:=sum(X[1:nSfunfun])*gs + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*gl;
annotation(smoothOrder=5);
end calc_g;
