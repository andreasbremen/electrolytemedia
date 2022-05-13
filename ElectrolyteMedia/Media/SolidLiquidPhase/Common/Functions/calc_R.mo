within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_R
  input Modelica.Media.Interfaces.Types.MassFraction[nSfunfun + nLfunfun] X;
  output SI.SpecificEntropy R;
protected
  SI.SpecificEntropy Rs = SolidFunctions.calc_R(X[1:nSfunfun]);
  SI.SpecificEntropy Rl = LiquidFunctions.calc_R(X[1+nSfunfun:nSfunfun+nLfunfun]);
algorithm
  R :=Rs + Rl;
  annotation(Inline=true,smoothOrder=5);
end calc_R;
