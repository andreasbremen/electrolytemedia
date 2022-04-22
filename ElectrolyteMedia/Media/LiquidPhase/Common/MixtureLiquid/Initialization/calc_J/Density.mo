within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid.Initialization.calc_J;
function Density
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Pressure p;
  input MassFraction[nL] X;
  output Real[1,1] J;
protected
  SI.DerDensityByPressure ddpT=Functions.calc_der_d_p(
      T,
      p,
      X);
algorithm
  J[1,1] :=ddpT;

end Density;
