within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.tautau;
function calc_Phi_Phi_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] Phi_Phi_ijtautau;
protected
  Real I = calc_I(X);
  Real[nLifun,nLifun] Phi_ijtautau= GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau.calc_Phi_ijtautau(T,p,X);
  Real[nLifun,nLifun] der_Phi_ijtautau= GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau.calc_der_E_theta_ijtautau(T,p,X);
algorithm

  Phi_Phi_ijtautau :=Phi_ijtautau + I*der_Phi_ijtautau;
end calc_Phi_Phi_ijtautau;
