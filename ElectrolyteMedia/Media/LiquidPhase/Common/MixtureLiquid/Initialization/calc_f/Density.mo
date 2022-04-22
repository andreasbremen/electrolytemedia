within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid.Initialization.calc_f;
function Density
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input SI.Temperature T;
  input SI.Pressure p;
  input MassFraction[nL] X;
  input SI.Density d;
  output Real[1] f;

// protected
//   SI.MoleFraction[nL] Yfull = x[1:nL]/sum(x[1:nL]);
//   SI.MassFraction[nL] Xfull = Functions.calc_X(Yfull,data,nL);
algorithm
  //isopotential and mass balance from LE_Tp
//   f[1:nL] :=LE_Tp(x[1:nL],c_i,T,x[nL + 1]);
  //Enthalpy balance
  f[1] :=(Functions.calc_d(
          T,
          p,
          X) - d);
end Density;
