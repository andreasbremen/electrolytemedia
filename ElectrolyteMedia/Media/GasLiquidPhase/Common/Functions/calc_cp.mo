within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_cp
  "calculates specific heat capacity cp of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Density dg;
  input SI.MassFraction[nGfunfun+nLfunfun] X;

  output SI.SpecificHeatCapacityAtConstantPressure cp;

protected
  SI.MassFraction[nGfunfun] Xg = X[1:nGfunfun]/sum(X[1:nGfunfun]);
  SI.Pressure p = GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nLfunfun] Xl = X[1+nGfunfun:nGfunfun+nLfunfun]/sum(X[1+nGfunfun:nGfunfun+nLfunfun]);
  SI.SpecificHeatCapacityAtConstantPressure cpg = GasFunctions.calc_cp(T,dg,Xg);
  SI.SpecificHeatCapacityAtConstantPressure cpl = LiquidFunctions.calc_cp(T,p,Xl);
algorithm

  cp:=sum(X[1:nGfunfun])*cpg + sum(X[1+nGfunfun:nGfunfun+nLfunfun])*cpl;

end calc_cp;
