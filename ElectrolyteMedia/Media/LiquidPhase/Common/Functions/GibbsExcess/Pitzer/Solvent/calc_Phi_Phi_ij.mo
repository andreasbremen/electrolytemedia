within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solvent;
function calc_Phi_Phi_ij
  "Calculates cation i - anion j interaction B_Phi_ij"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] Phi_Phi_ij;
protected
  Real I = calc_I(X);
  Real[nLifun,nLifun] Phi_ij= Solute.calc_Phi_ij(T,p,X);
  Real[nLifun,nLifun] der_Phi_ij= Solute.calc_der_Phi_ij(T,p,X);
algorithm

  Phi_Phi_ij :=Phi_ij + I*der_Phi_ij;
end calc_Phi_Phi_ij;
