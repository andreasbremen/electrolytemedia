within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Solvent;
function calc_ln_a
  "calculates decadic logarithm of water activity based on extended Debye Hueckel model"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real ln_a;
protected
  Real Phi = calc_Phi(T,p,X);
  Real[nLifun] mol_i = calc_mfromX(X);
  Real mol_s = sum(mol_i);
algorithm
  ln_a :=-Phi*IF97.MH2O*mol_s;

end calc_ln_a;
