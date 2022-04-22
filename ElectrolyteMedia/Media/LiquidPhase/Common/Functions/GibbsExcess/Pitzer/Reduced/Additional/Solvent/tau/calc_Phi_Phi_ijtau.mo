within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.tau;
function calc_Phi_Phi_ijtau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] Phi_Phi_ijtau;
protected
  Real I = calc_I(X);
  Real[nLifun,nLifun] Phi_ijtau= GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_Phi_ijtau(T,p,X);
  Real[nLifun,nLifun] der_Phi_ijtau= GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_der_E_theta_ijtau(T,p,X);
algorithm

  Phi_Phi_ijtau :=Phi_ijtau + I*der_Phi_ijtau;
end calc_Phi_Phi_ijtau;
