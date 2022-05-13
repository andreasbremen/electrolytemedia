within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_cv
  "calculates specific heat capacity cv of solid mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output SI.SpecificHeatCapacityAtConstantVolume cv;

protected
  SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
  SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
  SI.SpecificHeatCapacityAtConstantVolume cvs = SolidFunctions.calc_cp(T,Xs);
  SI.SpecificHeatCapacityAtConstantVolume cvl = LiquidFunctions.calc_cv(T,p,Xl);
algorithm

  cv:=sum(X[1:nSfunfun])*cvs + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*cvl;

end calc_cv;
