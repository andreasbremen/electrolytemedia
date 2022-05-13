within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_s
  "calculates specific entropy of gas, solid mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output SI.SpecificEntropy s;

protected
  SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
  SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
  SI.SpecificEntropy ss = SolidFunctions.calc_s(T,Xs);
  SI.SpecificEntropy sl = LiquidFunctions.calc_s(T,p,Xl);
algorithm

  s:=sum(X[1:nSfunfun])*ss + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*sl;

  annotation(smoothOrder=5);
end calc_s;
