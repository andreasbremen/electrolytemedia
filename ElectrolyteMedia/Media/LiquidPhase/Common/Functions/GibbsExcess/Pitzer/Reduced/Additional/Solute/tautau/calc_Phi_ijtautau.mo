within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_Phi_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] Phi_ijtautau;
protected
  Real[nLifun,nLifun] theta_ijtautau=calc_theta_ijtautau(
                                                   T);
  Real[nLifun,nLifun] E_theta_ijtautau=calc_E_theta_ijtautau(
      T,
      p,
      X);
algorithm
  Phi_ijtautau :=theta_ijtautau + E_theta_ijtautau;

end calc_Phi_ijtautau;
