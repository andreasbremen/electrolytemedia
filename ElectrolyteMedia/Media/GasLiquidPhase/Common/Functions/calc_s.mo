within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_s
  "calculates specific entropy of gas, solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Density dg;
  input SI.MassFraction[nGfunfun+nLfunfun] X;

  output SI.SpecificEntropy s;

protected
  SI.MassFraction[nGfunfun] Xg = X[1:nGfunfun]/sum(X[1:nGfunfun]);
  SI.Pressure p = GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nLfunfun] Xl = X[1+nGfunfun:nGfunfun+nLfunfun]/sum(X[1+nGfunfun:nGfunfun+nLfunfun]);
  SI.SpecificEntropy sg = GasFunctions.calc_s(T,dg,Xg);
  SI.SpecificEntropy sl = LiquidFunctions.calc_s(T,p,Xl);
algorithm

  s:=sum(X[1:nGfunfun])*sg + sum(X[1+nGfunfun:nGfunfun+nLfunfun])*sl;

  annotation(smoothOrder=5);
end calc_s;
