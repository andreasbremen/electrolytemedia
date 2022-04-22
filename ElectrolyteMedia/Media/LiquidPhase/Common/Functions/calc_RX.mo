within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_RX "Returns vector containing gas constants of all species"
  output SI.SpecificEntropy[nLfun] R;

algorithm
  for i in 1:nLfun-1 loop
    R[i] :=datafun[i].R;
  end for;
  R[nLfun] := IF97.RH2O;

end calc_RX;
