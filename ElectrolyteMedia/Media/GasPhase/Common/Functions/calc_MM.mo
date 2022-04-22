within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_MM "Returns molar mass of mixture"

  input SI.MassFraction[nGfun] X;
  output SI.MolarMass MM;

algorithm
    MM :=1/sum(X[j]/datafun[j].MM for j in 1:nGfun);

end calc_MM;
