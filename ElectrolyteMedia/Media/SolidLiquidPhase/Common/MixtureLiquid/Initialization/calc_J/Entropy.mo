within ElectrolyteMedia.Media.SolidLiquidPhase.Common.MixtureLiquid.Initialization.calc_J;
function Entropy
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input Temperature T;
  input SI.Pressure p;
  input MassFraction[ns+nL] X;
  output Real[1,1] J;

protected
  SI.SpecificHeatCapacity cp=Functions.calc_cp(T,p,X);
algorithm
  J[1,1] :=cp/T/Common.StdRefH2O.s_tr;

end Entropy;
