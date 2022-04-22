within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_der_Phi_ij "Unsymmetrical mixing, Pitzer p.123"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] der_Phi_ij;
protected
  Real[nLifun,nLifun] der_E_theta_ij=Solute.calc_der_E_theta_ij(
      T,
      p,
      X);
algorithm
  der_Phi_ij := der_E_theta_ij;

end calc_der_Phi_ij;
