within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid;
function calc_Xred
  input SI.MassFraction[nG+nL] Xfull;
  output SI.MassFraction[nX] Xred;

algorithm
  Xred :=transpose(lambda_mass)*Xfull/(sum(transpose(lambda_mass)*Xfull));

end calc_Xred;
