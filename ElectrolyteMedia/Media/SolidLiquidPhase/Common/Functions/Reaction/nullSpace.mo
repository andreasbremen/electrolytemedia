within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions.Reaction;
function nullSpace "Return the orthonormal nullspace of a matrix"
  extends Modelica.Icons.Function;

  input Real[:,:] nu;
  output Real[size(nu,2),size(nu,2)-size(nu,1)] nullSpace;
protected
  Integer nullity;
algorithm

  (nullSpace,nullity) :=Modelica.Math.Matrices.nullSpace(nu);

//   input Real A[nL, nY] "Input matrix";
//   input Integer nL;
//   input Integer nY;
//   output Real Z[size(A, 2), nY-nL] "Orthonormal nullspace of matrix A";
//   //output Integer nullity "Nullity, i.e., the dimension of the nullspace";
//
// protected
//   Real V[size(A, 2), size(A, 2)] "Right orthogonal matrix";
//   Real sigma[min(size(A, 1), size(A, 2))] "singular values";
//   Integer rank "rank of matrix A";
//   Real eps "tolerance for rank determination";
//   Integer n=min(size(A, 1), size(A, 2));
//   Integer i=n;
//   Integer nullity;
//
// algorithm
//   (sigma,,V) := Modelica.Math.Matrices.singularValues(A);
//   V := transpose(V);
//   // rank computation
//   eps := max(size(A, 1), size(A, 2))*max(sigma)*Modelica.Constants.eps;
//   rank := 0;
//   if n > 0 then
//     while i > 0 loop
//       if sigma[i] > eps then
//         rank := i;
//         i := 0;
//       end if;
//       i := i - 1;
//     end while;
//   end if;
//   Z := V[:, rank + 1:size(A, 2)];
//   // nullspace computation
//   nullity := size(A, 2) - rank;
//   // nullity
end nullSpace;
