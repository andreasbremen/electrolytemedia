within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_cv
  "calculates specific heat capacity cv of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Density dg;
  input SI.MassFraction[nGfunfun+nLfunfun] X;

  output SI.SpecificHeatCapacityAtConstantVolume cv;

protected
  SI.MassFraction[nGfunfun] Xg = X[1:nGfunfun]/sum(X[1:nGfunfun]);
  SI.Pressure p = GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nLfunfun] Xl = X[1+nGfunfun:nGfunfun+nLfunfun]/sum(X[1+nGfunfun:nGfunfun+nLfunfun]);
  SI.SpecificHeatCapacityAtConstantVolume cvg = GasFunctions.calc_cv(T,dg,Xg);
  SI.SpecificHeatCapacityAtConstantVolume cvl = LiquidFunctions.calc_cv(T,p,Xl);
algorithm

  cv:=sum(X[1:nGfunfun])*cvg + sum(X[1+nGfunfun:nGfunfun+nLfunfun])*cvl;

end calc_cv;
