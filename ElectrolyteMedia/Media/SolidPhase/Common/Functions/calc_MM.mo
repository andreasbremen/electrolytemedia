within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_MM "Mixture molar mass"

  input SI.MassFraction[nSfun] X;
  output SI.MolarMass MM;

algorithm
    MM :=1/(sum(X[j]/datafun[j].MM for j in 1:nSfun));

end calc_MM;
