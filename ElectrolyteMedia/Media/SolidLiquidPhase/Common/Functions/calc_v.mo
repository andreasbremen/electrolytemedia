within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_v
  "calculates specific volume of solid mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output SI.SpecificVolume v;

protected
  SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
  SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
  SI.SpecificVolume vs = SolidFunctions.calc_v(T,p,Xs);
  SI.SpecificVolume vl = LiquidFunctions.calc_v(T,p,Xl);
algorithm

  v:=sum(X[1:nSfunfun])*vs + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*vl;
  annotation(smoothOrder=5);
end calc_v;
