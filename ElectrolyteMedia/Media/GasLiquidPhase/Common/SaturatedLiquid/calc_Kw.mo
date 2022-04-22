within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid;
function calc_Kw
  "calculates equilibrium constant of water at temperature and pressure"
  input SI.Temperature T;
  input SI.Pressure p;
  output Real Kw;

protected
  SI.SpecificEnergy[nG+nL] gi;
  SI.MolarEnergy[nG+nL] gim;
  SI.MolarEnergy grw;
algorithm
  gi :=Functions.calc_gi(T, p);
  for i in 1:nG+nL loop
    gim[i] :=gi[i]*MMX[i];
  end for;
  grw :=gim[nG + nL] - gim[nG];
  Kw :=exp(-grw/Modelica.Constants.R/T);

end calc_Kw;
