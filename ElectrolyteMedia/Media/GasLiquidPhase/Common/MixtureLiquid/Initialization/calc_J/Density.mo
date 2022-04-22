within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid.Initialization.calc_J;
function Density
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Density dg;
  input MassFraction[nG+nL] X;
  output Real[1,1] J;
protected
  SI.DerDensityByPressure ddpT= Functions.calc_der_d_p(T,dg,X);
  SI.Density d=Functions.calc_d(T,dg,X);
algorithm
   J[1,1] :=ddpT;

end Density;
