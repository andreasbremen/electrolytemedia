within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_g
  "calculates specific Gibbs free energy of gas, solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Density dg;
  input SI.MassFraction[nGfunfun+nLfunfun] X;

  output SI.SpecificGibbsFreeEnergy g;

protected
  SI.MassFraction[nGfunfun] Xg = X[1:nGfunfun]/sum(X[1:nGfunfun]);
  SI.Pressure p = GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nLfunfun] Xl = X[1+nGfunfun:nGfunfun+nLfunfun]/sum(X[1+nGfunfun:nGfunfun+nLfunfun]);
  SI.SpecificGibbsFreeEnergy gg = GasFunctions.calc_g(T,dg,Xg);
  SI.SpecificGibbsFreeEnergy gl = LiquidFunctions.calc_g(T,p,Xl);
algorithm

  g:=sum(X[1:nGfunfun])*gg + sum(X[1+nGfunfun:nGfunfun+nLfunfun])*gl;

//   annotation(smoothOrder(normallyConstant=data)=2);
end calc_g;
