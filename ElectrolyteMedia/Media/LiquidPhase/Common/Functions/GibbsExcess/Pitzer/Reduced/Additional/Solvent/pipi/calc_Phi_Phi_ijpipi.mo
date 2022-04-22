within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.pipi;
function calc_Phi_Phi_ijpipi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] Phi_Phi_ijpipi;
protected
  Real I = calc_I(X);
  Real[nLifun,nLifun] Phi_ijpipi= GibbsExcess.Pitzer.Reduced.Additional.Solute.pipi.calc_Phi_ijpipi(T,p,X);
  Real[nLifun,nLifun] der_Phi_ijpipi= GibbsExcess.Pitzer.Reduced.Additional.Solute.pipi.calc_der_E_theta_ijpipi(T,p,X);
algorithm

  Phi_Phi_ijpipi :=Phi_ijpipi + I*der_Phi_ijpipi;
end calc_Phi_Phi_ijpipi;
