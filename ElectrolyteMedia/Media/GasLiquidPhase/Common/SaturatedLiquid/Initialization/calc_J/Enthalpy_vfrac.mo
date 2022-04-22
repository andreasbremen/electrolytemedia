within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid.Initialization.calc_J;
function Enthalpy_vfrac
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Density dg;
  input MassFraction[nG+nL] X;
  output Real[1,1] J;
protected
  SI.MassFraction[nG] Xg = X[1:nG]/sum(X[1:nG]);
  SI.Pressure p = Functions.GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nL] Xl = X[1+nG:nG+nL]/sum(X[1+nG:nG+nL]);
  SI.SpecificEnthalpy hg = Functions.GasFunctions.calc_h(T,dg,Xg);
  SI.SpecificEnthalpy hl = Functions.LiquidFunctions.calc_h(T,p,Xl);
algorithm
  J[1,1] :=(hg - hl)/Common.StdRefH2O.h_tr;

end Enthalpy_vfrac;
