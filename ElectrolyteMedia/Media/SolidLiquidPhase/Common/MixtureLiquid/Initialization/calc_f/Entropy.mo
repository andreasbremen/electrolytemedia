within ElectrolyteMedia.Media.SolidLiquidPhase.Common.MixtureLiquid.Initialization.calc_f;
function Entropy
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.SpecificEntropy s;
  input MassFraction[ns+nL] X;
  output Real[1] f;

algorithm
  f[1] :=(Functions.calc_s(T,p,X) - s)/Common.StdRefH2O.s_tr;
end Entropy;
