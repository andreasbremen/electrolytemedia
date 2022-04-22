within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_MMX "Returns vector containing molar masses of all species"
  output SI.MolarMass[nLfun] MMX;

algorithm
  MMX[1:nLifun] :=datafun[1:nLifun].MM;
  MMX[nLfun] := IF97.MH2O;
  annotation(smoothOrder=5);
end calc_MMX;
