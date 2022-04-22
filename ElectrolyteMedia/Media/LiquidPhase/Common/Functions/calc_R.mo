within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_R "Mass based gas constant"
  input SI.MassFraction[nLfun] X;
  output SI.SpecificEntropy R;

protected
  SI.SpecificEntropy[nLfun] Ri;
algorithm
  for i in 1:nLfun-1 loop
    Ri[i] :=datafun[i].R;
  end for;
  Ri[nLfun] := IF97.RH2O;

  R :=X*Ri;

  annotation(smoothOrder=5);
end calc_R;
