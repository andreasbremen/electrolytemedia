within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pitau;
function calc_Phi_ijpitau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] Phi_ijpitau;
protected
  Real[nLifun,nLifun] E_theta_ijpitau=calc_E_theta_ijpitau(
      T,
      p,
      X);
algorithm
  Phi_ijpitau := E_theta_ijpitau;

end calc_Phi_ijpitau;
