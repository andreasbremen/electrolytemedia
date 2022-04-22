within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_MM
  "calculates molar mass of gas, solute and solvent mixture at T and p"
  input SI.MassFraction[nGfunfun+nLfunfun] X;

  output SI.MolarMass MM;

protected
  SI.MassFraction[nGfunfun] Xg = X[1:nGfunfun]/sum(X[1:nGfunfun]);
  SI.MassFraction[nLfunfun] Xl = X[1+nGfunfun:nGfunfun+nLfunfun]/sum(X[1+nGfunfun:nGfunfun+nLfunfun]);
  SI.MolarMass MMg = GasFunctions.calc_MM(Xg);
  SI.MolarMass MMl = LiquidFunctions.calc_MM(Xl);
algorithm
    MM :=1/(sum(X[1:nGfunfun])/MMg + sum(X[1+nGfunfun:nGfunfun+nLfunfun])/MMl);

end calc_MM;
