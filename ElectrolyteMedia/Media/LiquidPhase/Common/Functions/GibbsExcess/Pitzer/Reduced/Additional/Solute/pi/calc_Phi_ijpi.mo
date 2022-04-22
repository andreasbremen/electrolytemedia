within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pi;
function calc_Phi_ijpi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] Phi_ijpi;
protected
  Real[nLifun,nLifun] E_theta_ijpi=calc_E_theta_ijpi(
      T,
      p,
      X);
algorithm
  Phi_ijpi := E_theta_ijpi;

end calc_Phi_ijpi;
