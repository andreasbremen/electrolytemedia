within ElectrolyteMedia.Media.Common.Reaction;
function calc_lambda_id
  "calculates lambda with nX*nX identity matrix"
  input Real[:,:] lambda;

  output Real[size(lambda,1),size(lambda,2)] lambda_id;
  output Integer[size(lambda,1),size(lambda,1)] P_id;
protected
  Integer nF = size(lambda,1);
  Integer nX = size(lambda,2);
  Real[nF,nX] LU;
  Integer[nX] pivots_LAPACK;
  Real[nF,nX] L;
  Integer[nF,nF] P;
  Real[nX,nX] Lsquare;
  Real[nX,nX] LsquareT;
  Real[nX,nX] LU_;
  Integer[nX] pivots_LAPACK_;
  Real[nX,nX] U_;

algorithm
  (LU,pivots_LAPACK) :=Modelica.Math.Matrices.LU(lambda);
  L :=Media.Common.Reaction.calc_L(LU);
  P :=Media.Common.Reaction.calc_P(pivots_LAPACK,nF);
  Lsquare :=L[1:nX, :];
  LsquareT :=transpose(Lsquare);
  (LU_,pivots_LAPACK_) :=Modelica.Math.Matrices.LU(LsquareT);
  U_ :=Media.Common.Reaction.calc_U(LU_);

  lambda_id :=L*Modelica.Math.Matrices.inv(transpose(U_));
  P_id :=P;
end calc_lambda_id;
