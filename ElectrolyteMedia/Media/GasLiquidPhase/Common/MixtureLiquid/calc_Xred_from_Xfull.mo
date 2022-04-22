within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid;
function calc_Xred_from_Xfull
  "Calculates reduced mass fraction vector from full mass fraction vector"
  input SI.MassFraction[nF] Xfull;
  output SI.MassFraction[nX] Xred;

protected
  SI.Mass[nX] mred = transpose(lambda_mass)*Xfull;
  SI.Mass mtotal = sum(mred);
algorithm
  Xred :=mred/mtotal;
//   Xred :=if nR > 0 then transpose(lambda_mass)*Xfull/(sum(transpose(lambda_mass)*Xfull)) else Xfull;

end calc_Xred_from_Xfull;
