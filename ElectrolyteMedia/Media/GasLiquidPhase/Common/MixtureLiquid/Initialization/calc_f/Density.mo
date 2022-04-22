within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid.Initialization.calc_f;
function Density
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Density dg;
  input MassFraction[nG+nL] X;
  input SI.Density d;
  output Real[1] f;

algorithm
//    f[1] := (log(Functions.calc_d(T,p,X,datag,datal,nG,nL)) - log(d));
   f[1] :=(Functions.calc_d(
    T,
    dg,
    X) - d);
end Density;
