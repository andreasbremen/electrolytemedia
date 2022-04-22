within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid.Initialization.calc_f;
function Enthalpy
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.SpecificEnthalpy h;
  input MassFraction[nL] X;
  output Real[1] f;

// protected
//   SI.MoleFraction[nL] Yfull = x[1:nL]/sum(x[1:nL]);
//   SI.MassFraction[nL] Xfull = Functions.calc_X(Yfull,data,nL);
algorithm
  //Enthalpy balance
  f[1] :=(Functions.calc_h(
          T,
          p,
          X) - h)/Common.StdRefH2O.h_tr;                      //298.15-x[nL+1];//(Functions.calc_h(x[nL+1],p,Xfull,data,nL) - h)*1e-8;
end Enthalpy;
