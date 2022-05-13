within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_MM
  "calculates molar mass of gas, solid mixture at T and p"
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output SI.MolarMass MM;

protected
  SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
  SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
  SI.MolarMass MMs = SolidFunctions.calc_MM(Xs);
  SI.MolarMass MMl = LiquidFunctions.calc_MM(Xl);
algorithm
    MM :=1/(sum(X[1:nSfunfun])/MMs + sum(X[1+nSfunfun:nSfunfun+nLfunfun])/MMl);

end calc_MM;
