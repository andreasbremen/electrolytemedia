within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_mfromX
  "Calculates molality vector from mass fraction vector of liquid species"
  input SI.MassFraction[:] X;
  output SI.MoleFraction[nLfun-1] m;
protected
  SI.MoleFraction[nLfun] Y = calc_Y(X);

algorithm

  m :=calc_mfromY(Y);

  annotation(Inline=true,smoothOrder=5);
end calc_mfromX;
