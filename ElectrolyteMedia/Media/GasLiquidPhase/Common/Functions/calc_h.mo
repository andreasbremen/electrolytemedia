within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_h
  "calculates specific enthalpy of gas, solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Density dg "Gas phase density";
  input SI.MassFraction[nGfunfun+nLfunfun] X;

  output SI.SpecificEnthalpy h;

protected
  SI.MassFraction[nGfunfun] Xg = X[1:nGfunfun]/sum(X[1:nGfunfun]);
  SI.Pressure p = GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nLfunfun] Xl = X[1+nGfunfun:nGfunfun+nLfunfun]/sum(X[1+nGfunfun:nGfunfun+nLfunfun]);
  SI.SpecificEnthalpy hg = GasFunctions.calc_h(T,dg,Xg);
  SI.SpecificEnthalpy hl = LiquidFunctions.calc_h(T,p,Xl);
algorithm

  h:=sum(X[1:nGfunfun])*hg + sum(X[1+nGfunfun:nGfunfun+nLfunfun])*hl;

  annotation(smoothOrder=5);
end calc_h;
