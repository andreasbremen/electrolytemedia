within ElectrolyteMedia.Media.Common.Reaction;
function calc_P
  "calculates permuation matrix of pivot vector from LAPACK"
  input Integer[:] pivots_LAPACK_;
  input Integer nF;
  output Integer [nF,nF] P;
protected
  Integer nR = size(pivots_LAPACK_,1);
  Integer[nF] pivots;
  Integer[nF] pivots_LAPACK;
  Integer temp;
algorithm
  // convert pivots_LAPACK_ to pivot vector
  pivots :={i for i in 1:nF};
  pivots_LAPACK[1:nR] :=pivots_LAPACK_;
  pivots_LAPACK[nR+1:nF] :={i for i in nR + 1:nF};
  for i in 1:nF-1 loop
    temp :=pivots[pivots_LAPACK[i]];
    pivots[pivots_LAPACK[i]] :=pivots[i];
    pivots[i] :=temp;
  end for;

  // compute permutation matrix P
  for i in 1:nF loop
    P[pivots[i],i] :=1;
  end for;

end calc_P;
