within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_u
  "calculates specific internal energy of gas, solid mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output SI.SpecificInternalEnergy u;

protected
  SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
  SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
  SI.SpecificInternalEnergy us = SolidFunctions.calc_u(T,p,Xs);
  SI.SpecificInternalEnergy ul = LiquidFunctions.calc_u(T,p,Xl);
algorithm

  u:=sum(X[1:nSfunfun])*us + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*ul;

  annotation(Inline=true,smoothOrder=5);
end calc_u;
