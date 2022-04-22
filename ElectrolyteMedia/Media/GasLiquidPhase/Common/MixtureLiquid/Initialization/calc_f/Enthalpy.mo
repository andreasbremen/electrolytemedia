within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid.Initialization.calc_f;
function Enthalpy
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Density dg;
  input SI.SpecificEnthalpy h;
  input MassFraction[nG+nL] X;
  output Real[1] f;

algorithm
  //Enthalpy balance
  f[1] :=(Functions.calc_h(T,dg,X) - h)/Common.StdRefH2O.h_tr;
end Enthalpy;
