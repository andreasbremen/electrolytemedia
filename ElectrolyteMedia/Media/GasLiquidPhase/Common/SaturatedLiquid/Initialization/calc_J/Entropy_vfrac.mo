within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid.Initialization.calc_J;
function Entropy_vfrac
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Density dg;
  input MassFraction[nG+nL] X;
  output Real[1,1] J;
protected
  SI.MassFraction[nG] Xg = X[1:nG]/sum(X[1:nG]);
  SI.Pressure p = Functions.GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nL] Xl = X[1+nG:nG+nL]/sum(X[1+nG:nG+nL]);
  SI.SpecificEntropy sg = Functions.GasFunctions.calc_s(T,dg,Xg);
  SI.SpecificEntropy sl = Functions.LiquidFunctions.calc_s(T,p,Xl);
algorithm
  J[1,1] :=(sg - sl)/Common.StdRefH2O.s_tr;

end Entropy_vfrac;
