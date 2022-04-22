within ElectrolyteMedia.Media.Common.Reaction;
function calc_L "calculates lower L of LU decomposition"
  input Real[:,:] LU;
  output Real[size(LU,1),size(LU,2)] L;
protected
  Integer nF = size(LU,1);
  Integer nR = size(LU,2);
algorithm
  //compute lower L
  for i in 1:nF loop
    for j in 1:nR loop
      if i == j then
        L[i,j] :=1;
      elseif j < i then
        L[i,j] :=LU[i, j];
      end if;
    end for;
  end for;
end calc_L;
