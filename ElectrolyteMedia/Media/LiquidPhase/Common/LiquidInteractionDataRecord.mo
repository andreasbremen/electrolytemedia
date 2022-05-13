within ElectrolyteMedia.Media.LiquidPhase.Common;
record LiquidInteractionDataRecord "Coefficient data record for interaction properties between aqueous solute species"

  parameter Integer nLi = 1;

  //Bromley parameter
  Real[nLi,nLi] Bromley_ij = zeros(nLi,nLi) "Bromley parameter";

  //Pitzer parameters
   Real [nLi,nLi,8] beta0_ij = zeros(nLi,nLi,8);
   Real [nLi,nLi,8] beta1_ij = zeros(nLi,nLi,8);
   Real [nLi,nLi,6] beta2_ij = zeros(nLi,nLi,6);
   Real [nLi,nLi,8] c_ij = zeros(nLi,nLi,8);
   Real [nLi,nLi,6] theta_ij = zeros(nLi,nLi,6);
   Real [nLi,nLi,nLi,6] psi_ijk = zeros(nLi,nLi,nLi,6);
   Real [nLi,nLi,8] lambda_ij = zeros(nLi,nLi,8);
   Real [nLi,nLi,nLi,8] zeta_ijk = zeros(nLi,nLi,nLi,8);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));

end LiquidInteractionDataRecord;
