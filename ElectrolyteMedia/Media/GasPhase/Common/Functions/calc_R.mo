within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_R "Mass based gas constant"
  input SI.MassFraction[nGfun] X;

  output SI.SpecificEntropy R;

algorithm
    R :=datafun[:].R*X;
annotation(smoothOrder=5);
end calc_R;
