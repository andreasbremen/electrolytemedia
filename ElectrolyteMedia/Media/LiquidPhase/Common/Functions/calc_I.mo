within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_I
  "Calculates ionic strength in molality basis"

  input Real[nLfun] X "mass fraction";
  output Real I "Ionic strength";

protected
  Real[nLfun - 1] mol=calc_mfromX(X);

algorithm
  I := 0.5 * sum(mol[i] * datafun[i].z^2 for i in 1:nLfun-1);
  annotation(Inline=true,smoothOrder = 5);
end calc_I;
