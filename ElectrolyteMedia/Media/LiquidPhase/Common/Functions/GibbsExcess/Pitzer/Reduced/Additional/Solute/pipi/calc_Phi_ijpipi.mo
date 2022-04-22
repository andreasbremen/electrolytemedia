within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pipi;
function calc_Phi_ijpipi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] Phi_ijpipi;
protected
  Real[nLifun,nLifun] E_theta_ijpipi=calc_E_theta_ijpipi(
      T,
      p,
      X);
algorithm
  Phi_ijpipi := E_theta_ijpipi;

end calc_Phi_ijpipi;
