within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_v
  "calculates specific volume of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Density dg;
  input SI.MassFraction[nGfunfun+nLfunfun] X;

  output SI.SpecificVolume v;

protected
  SI.MassFraction[nGfunfun] Xg = X[1:nGfunfun]/sum(X[1:nGfunfun]);
  SI.Pressure p = GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nLfunfun] Xl = X[1+nGfunfun:nGfunfun+nLfunfun]/sum(X[1+nGfunfun:nGfunfun+nLfunfun]);
  SI.SpecificVolume vg = 1/dg;//GasFunctions.calc_v(T,p,Xg);
  SI.SpecificVolume vl = LiquidFunctions.calc_v(T,p,Xl);
algorithm

  v:=sum(X[1:nGfunfun])*vg + sum(X[1+nGfunfun:nGfunfun+nLfunfun])*vl;
  annotation(smoothOrder=5);
end calc_v;
