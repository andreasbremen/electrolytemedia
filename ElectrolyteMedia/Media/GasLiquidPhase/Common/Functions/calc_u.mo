within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_u
  "calculates specific internal energy of gas, solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Density dg;
  input SI.MassFraction[nGfunfun+nLfunfun] X;

  output SI.SpecificInternalEnergy u;

protected
  SI.MassFraction[nGfunfun] Xg = X[1:nGfunfun]/sum(X[1:nGfunfun]);
  SI.Pressure p = GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nLfunfun] Xl = X[1+nGfunfun:nGfunfun+nLfunfun]/sum(X[1+nGfunfun:nGfunfun+nLfunfun]);
  SI.SpecificInternalEnergy ug = GasFunctions.calc_u(T,dg,Xg);
  SI.SpecificInternalEnergy ul = LiquidFunctions.calc_u(T,p,Xl);
algorithm

  u:=sum(X[1:nGfunfun])*ug + sum(X[1+nGfunfun:nGfunfun+nLfunfun])*ul;

  annotation(smoothOrder=5);
end calc_u;
