within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_d
  "calculates density of gas, solid mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output SI.Density d;
protected
  SI.SpecificVolume v;
algorithm

  v :=calc_v(T,p,X);
  d :=1/v;

  annotation(Inline=true,smoothOrder=5);
end calc_d;
