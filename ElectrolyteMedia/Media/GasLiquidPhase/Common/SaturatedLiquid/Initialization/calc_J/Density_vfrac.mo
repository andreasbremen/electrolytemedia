within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid.Initialization.calc_J;
function Density_vfrac
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Density dg;
  input MassFraction[nG+nL] X;
  output Real[1,1] J;
protected
  SI.MassFraction[nG] Xg = X[1:nG]/sum(X[1:nG]);
  SI.Pressure p = Functions.GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nL] Xl = X[1+nG:nG+nL]/sum(X[1+nG:nG+nL]);
  SI.SpecificVolume vg = 1/dg;
  SI.SpecificVolume vl = Functions.LiquidFunctions.calc_v(T,p,Xl);
  SI.SpecificVolume v = Functions.calc_v(T,dg,X);
algorithm
  J[1,1] :=-1/v^2*(vg - vl);

end Density_vfrac;
