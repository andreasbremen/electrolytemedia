within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_Phi_ijtau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] Phi_ijtau;
protected
  Real[nLifun,nLifun] theta_ijtau=calc_theta_ijtau(T);
  Real[nLifun,nLifun] E_theta_ijtau=calc_E_theta_ijtau(
      T,
      p,
      X);
algorithm
  Phi_ijtau :=theta_ijtau + E_theta_ijtau;

end calc_Phi_ijtau;
