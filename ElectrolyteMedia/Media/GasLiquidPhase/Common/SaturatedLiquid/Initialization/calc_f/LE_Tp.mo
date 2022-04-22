within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid.Initialization.calc_f;
function LE_Tp
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input Real[nL] x;
  input Real[nX_L] Xredl;
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nL] f;

protected
  SI.MolarEnergy gi[nL];
  SI.MolarEnergy gr[nR_L];
  Real[nR_L] logK;
  Real[nL] gamma;
  Real[nL] m;
  Real[nL] loga;
  SI.MoleFraction[nL] Y;
  SI.MassFraction[nL] X;
  SI.Mass[nL] mass;
algorithm

  Y :=x/sum(x);
  m[1:nL-1] := Functions.LiquidFunctions.calc_mfromY(Y);
  m[nL] := Y[nL]/MH2O;
  X :=Functions.LiquidFunctions.calc_X(Y);
  gamma := Functions.LiquidFunctions.calc_gamma(T,p,X);
  gi := Functions.LiquidFunctions.calc_gi(T, p) .* Functions.LiquidFunctions.calc_MMX();

  gr :=nu_L*gi;
  logK :=-gr/Modelica.Constants.R/T;

  loga :=log(gamma .* m);
  mass :=x .* Functions.LiquidFunctions.calc_MMX();

  f[1:nR_L] :=nu_L*loga - logK;
  f[nR_L+1:nL] :=transpose(lambda_mass_L)*mass - Xredl;

//   for j in 1:nR_L loop
//     gr[j] := sum({nu_L[j, i]*gi[i] for i in 1:nL});
//     K[j] :=exp(-gr[j]/Modelica.Constants.R/T);
//   end for;
//   //isopotential
//   for r in 1:nR_L loop
//     f[r] := log(product({(gamma[i_]*m[i_])^nu_L[r, i_] for i_ in 1:nL})) - log(K[r]);
//   end for;
//   //reduced mass balance
//   for j in 1:nL - nR_L loop
//     f[j + nR_L] := x[1:nL]*lambda_L[1:nL, j] - c_i[j];
//   end for;
end LE_Tp;
