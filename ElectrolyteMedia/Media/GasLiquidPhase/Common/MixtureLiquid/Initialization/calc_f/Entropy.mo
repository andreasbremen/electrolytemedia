within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid.Initialization.calc_f;
function Entropy
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Density dg;
  input SI.SpecificEntropy s;
  input MassFraction[nG+nL] X;
  output Real[1] f;

algorithm
  f[1] :=(Functions.calc_s(
    T,
    dg,
    X) - s)/Common.StdRefH2O.s_tr;
end Entropy;
