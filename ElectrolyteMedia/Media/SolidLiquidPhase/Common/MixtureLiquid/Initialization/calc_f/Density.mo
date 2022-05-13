within ElectrolyteMedia.Media.SolidLiquidPhase.Common.MixtureLiquid.Initialization.calc_f;
function Density
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Pressure p;
  input MassFraction[ns+nL] X;
  input SI.Density d;
  output Real[1] f;

algorithm
   f[1] :=(Functions.calc_d(T,p,X) - d);
end Density;
