within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_R
  input Modelica.Media.Interfaces.Types.MassFraction[nGfunfun + nLfunfun] X;
  output SI.SpecificEntropy R;
protected
  SI.SpecificEntropy Rg = GasFunctions.calc_R(X[1:nGfunfun]);
  SI.SpecificEntropy Rl = LiquidFunctions.calc_R(X[1+nGfunfun:nGfunfun+nLfunfun]);
algorithm
  R :=Rg + Rl;
end calc_R;
