within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.tau;
function calc_ln_atau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real ln_atau;
protected
  Real Phitau = calc_Phitau(T,p,X);
  Real[nLifun] mol_i = calc_mfromX(X);
  Real mol_s = sum(mol_i);
algorithm
  ln_atau :=-Phitau*IF97.MH2O*mol_s;

end calc_ln_atau;
