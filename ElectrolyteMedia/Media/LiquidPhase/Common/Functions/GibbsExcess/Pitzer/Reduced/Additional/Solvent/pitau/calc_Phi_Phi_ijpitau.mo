within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.pitau;
function calc_Phi_Phi_ijpitau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] Phi_Phi_ijpitau;
protected
  Real I = calc_I(X);
  Real[nLifun,nLifun] Phi_ijpitau= GibbsExcess.Pitzer.Reduced.Additional.Solute.pitau.calc_Phi_ijpitau(T,p,X);
  Real[nLifun,nLifun] der_Phi_ijpitau= GibbsExcess.Pitzer.Reduced.Additional.Solute.pitau.calc_der_E_theta_ijpitau(T,p,X);
algorithm

  Phi_Phi_ijpitau :=Phi_ijpitau + I*der_Phi_ijpitau;
end calc_Phi_Phi_ijpitau;
