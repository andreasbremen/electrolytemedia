within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_MM "Mixture molar mass"

  input SI.MassFraction[nLfun] X;
  output SI.MolarMass MM;

algorithm
    MM :=1/(sum(X[j]/datafun[j].MM for j in 1:nLfun-1)+X[nLfun]/IF97.MH2O);

end calc_MM;
