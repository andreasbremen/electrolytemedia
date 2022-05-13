within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_h
  "calculates specific enthalpy of gas, solid mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output SI.SpecificEnthalpy h;

protected
  SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
  SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
  SI.SpecificEnthalpy hs = SolidFunctions.calc_h(T,p,Xs);
  SI.SpecificEnthalpy hl = LiquidFunctions.calc_h(T,p,Xl);
algorithm

  h:=sum(X[1:nSfunfun])*hs + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*hl;

  annotation(Inline=true,smoothOrder=5);
end calc_h;
