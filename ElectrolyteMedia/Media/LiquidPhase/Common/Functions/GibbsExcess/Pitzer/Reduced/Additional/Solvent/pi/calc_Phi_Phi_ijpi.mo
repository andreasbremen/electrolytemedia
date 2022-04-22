within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.pi;
function calc_Phi_Phi_ijpi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] Phi_Phi_ijpi;
protected
  Real I = calc_I(X);
  Real[nLifun,nLifun] Phi_ijpi= GibbsExcess.Pitzer.Reduced.Additional.Solute.pi.calc_Phi_ijpi(T,p,X);
  Real[nLifun,nLifun] der_Phi_ijpi= GibbsExcess.Pitzer.Reduced.Additional.Solute.pi.calc_der_E_theta_ijpi(T,p,X);
algorithm

  Phi_Phi_ijpi :=Phi_ijpi + I*der_Phi_ijpi;
end calc_Phi_Phi_ijpi;
