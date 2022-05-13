within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_cp
  "calculates specific heat capacity cp of solid mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output SI.SpecificHeatCapacityAtConstantPressure cp;

protected
  SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
  SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
  SI.SpecificHeatCapacityAtConstantPressure cps = SolidFunctions.calc_cp(T,Xs);
  SI.SpecificHeatCapacityAtConstantPressure cpl = LiquidFunctions.calc_cp(T,p,Xl);
algorithm

  cp:=sum(X[1:nSfunfun])*cps + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*cpl;

end calc_cp;
