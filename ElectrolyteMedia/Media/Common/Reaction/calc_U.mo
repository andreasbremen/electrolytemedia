within ElectrolyteMedia.Media.Common.Reaction;
function calc_U "calculates upper U of LU decomposition"
  input Real[:,:] LU;
  output Real[size(LU,2),size(LU,2)] U;
protected
  Integer nF = size(LU,1);
  Integer nR = size(LU,2);
algorithm
  //compute upper U
   for i in 1:nR loop
    for j in 1:nR loop
      if i < j+1 then
        U[i,j] :=LU[i, j];
      end if;
    end for;
  end for;
end calc_U;
