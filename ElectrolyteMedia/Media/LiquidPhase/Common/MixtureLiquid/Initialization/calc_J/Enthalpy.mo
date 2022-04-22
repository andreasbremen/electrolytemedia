within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid.Initialization.calc_J;
function Enthalpy
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Pressure p;
  input MassFraction[nL] X;
  output Real[1,1] J;
protected
  SI.SpecificHeatCapacity cp=
      Functions.calc_cp(
            T,
            p,
            X);
algorithm
  J[1,1] :=cp/Common.StdRefH2O.h_tr;

end Enthalpy;
