within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid.Initialization.calc_f;
function LE_Tp
  "Liquid dissociation equilibrium with constraints on solute molalities"
  input Real[nL] x;
  input Real[nX] Xred;
  input SI.Temperature T;
  input SI.Pressure p;
  output Real[nL] f;

protected
  SI.MolarEnergy gi[nL];
  SI.MolarEnergy gr[nR];
  Real[nR] logK;
  Real[nL] gamma_i;
  Real prod;
  Real[nL] m;
  Real[nL] loga;
  SI.MoleFraction[nL] Y;
  SI.MassFraction[nL] X;
  SI.Mass[nL] mass;
algorithm

  Y :=x/sum(x);

  m[1:nL-1] := Functions.calc_mfromY(Y);
  m[nL] := Y[nL]/MH2O;

  X :=Functions.calc_X(Y);

  gamma_i := Functions.calc_gamma(T,p,X);

  gi := Functions.calc_gi(T, p) .* MMX;

  gr :=nu*gi;
  logK :=-gr/Modelica.Constants.R/T;

  loga :=log(gamma_i .* m);
  mass :=x .* MMX;

  //isopotential
  f[1:nR] :=nu*loga - logK;

  //reduced mass balance
  f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
end LE_Tp;
