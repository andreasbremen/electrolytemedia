within ElectrolyteMedia.Media.Common.Reaction;
function calc_nu_id
  "Calculates a stoichiometry matrix with pivoting so that nuout =(I|nu_var) with permutation matrix Pout for reordering species elements of nu_ordered = nuout*Pout"
  input Real[:,:] nu;

  output Real[size(nu,1),size(nu,2)] nuout;
  output Integer[size(nu,2),size(nu,2)] Pout;
protected
  parameter Integer nR = size(nu,1);
  parameter Integer nF = size(nu,2);
  Real[nF,nR] nuT = transpose(nu);
  Real[nR,nF] nuvar;
  Real[nF,nR] LU;
  Integer[nR] pivots_LAPACK;
  Integer info;
  Integer[nF] pivots;
  Real[nF,nR] L;
  Real[nR,nR] U;
  Integer[nF,nF] P;
  Real[nR,nR] Lsquare;
  Real[nR,nR] LsquareT;
  Real[nR,nR] LU_;
  Integer[nR] pivots_LAPACK_;
  Real[nR,nR] U_;
  Integer[nF,nF] P_;
  Real[nF,nR] nuT_;
algorithm

  (LU,pivots_LAPACK,info) :=Modelica.Math.Matrices.LU(nuT);

  // compute lower L from LU
  L :=calc_L(LU);

  // compute permutation matrix P from pivots_LAPACK
  P :=calc_P(pivots_LAPACK[1:nR],nF);

  // compute identity matrix of square of L such that L_ = (I|L') with L = P_*L_*U_

  Lsquare :=L[1:nR, :];
  LsquareT :=transpose(Lsquare);
  (LU_,pivots_LAPACK_) :=Modelica.Math.Matrices.LU(LsquareT);
  U_ :=calc_U(LU_);
  P_ :=identity(nF);
  P_[1:nR,1:nR] :=calc_P(pivots_LAPACK_, nR);

  nuT_ :=L*Modelica.Math.Matrices.inv(transpose(U_));

  nuout :=transpose(nuT_);
  Pout :=P*transpose(P_);

end calc_nu_id;
