within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_Phil "calculates liquid volume fraction"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;
  output SI.VolumeFraction Phil;
protected
  SI.MassFraction lfrac = sum(X[1+nSfunfun:nSfunfun+nLfunfun]);
  SI.MassFraction[nSfunfun] Xs = X[1:nSfunfun]/sum(X[1:nSfunfun]);
  SI.MassFraction[nLfunfun] Xl = X[1+nSfunfun:nSfunfun+nLfunfun]/sum(X[1+nSfunfun:nSfunfun+nLfunfun]);
  SI.SpecificVolume vl = LiquidFunctions.calc_v(T,p,Xl);
  SI.SpecificVolume vs = SolidFunctions.calc_v(T,p,Xs);
algorithm
  Phil :=lfrac*vl/(lfrac*vl + (1 - lfrac)*vs);

  annotation(Inline=true,smoothOrder=5);
end calc_Phil;
