within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_Phi_ij "Unsymmetrical mixing, Pitzer p.123"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] Phi_ij;
protected
  Real [nLifun,nLifun] theta_ij = calc_theta_ij(T);
  Real[nLifun,nLifun] E_theta_ij=Solute.calc_E_theta_ij(
      T,
      p,
      X);
algorithm
  Phi_ij :=theta_ij + E_theta_ij;

end calc_Phi_ij;
