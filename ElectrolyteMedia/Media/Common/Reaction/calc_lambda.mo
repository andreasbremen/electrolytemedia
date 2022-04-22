within ElectrolyteMedia.Media.Common.Reaction;
function calc_lambda
  "Calculates null vectors of stoichiometry matrix nu with positive elements only"
  input Real[:,:] nu;
  // input Real[nR,nF] nu;
  input Integer nR;
  input Integer nF;

  output Real[nF,nF-nR] lambda;

protected
  Integer[nF] sort;
  Real[nF,nR+nF] total;
  Real[1,nR+nF] total_ = -1*ones(1,nF+nR);
  Real[1,nR+nF] total__;
  Real[nF] iDummy;
  Integer temp;

  Real[nF-nR,nF] lambdaT;

  Real[nF] iDummy1;
  Integer[nF] sort1;
  Integer[nF,nF] resort1;
  Boolean[nF] resortindex = fill(false,nF);
  Integer index;
algorithm

  if nR > 0 then

  // total = (transpose(nu)|identity)
  total :=cat(2,transpose(nu),identity(nF));

  // Gauss Jordan algorithm with total
  for i in 1:nR loop
    if not total[i,i] <> 0 then
      (iDummy[i:nF],sort[i:nF]):=Modelica.Math.Vectors.sort(abs(total[i:nF, i]),false);
      for j in i:nF loop
        sort[j] :=sort[j] + (i-1);
      end for;
      total[i:nF,:] :=total[sort[i:nF], :];
    end if;
    total[i,:] :=total[i, :] / total[i, i];
    for j in 1:nF loop
      if i <> j then
        total[j,:] :=total[j, :] - total[j, i]*total[i, :]/total[i, i];
      end if;
    end for;
  end for;

  // read transpose(lambda) = lambdaT from total
  lambdaT :=total[nR + 1:nF, nR + 1:nR + nF];

  // column changes (saved in resort1) and row operations to get rid of negativ elements in 1:nF-nR columns of lambdaT
  for i in 1:nF-nR loop
    if min(lambdaT) < 0 then
      if not lambdaT[i,i] <> 0 then
        (iDummy1[i:nF],sort1[i:nF]) :=Modelica.Math.Vectors.sort(abs(lambdaT[i,i:nF]), false);
        for j in i:nF loop
          sort1[j] :=sort1[j] + (i - 1);
        end for;
        lambdaT[:,i:nF] :=lambdaT[:, sort1[i:nF]];
        resortindex[i] :=true;
        (iDummy1[i:nF],resort1[i:nF,i]) :=Modelica.Math.Vectors.sort(sort1[i:nF], true);
        for j in i:nF loop
          resort1[j,i] :=resort1[j,i] + (i-1);
        end for;
      end if;
       lambdaT[i,:] :=lambdaT[i, :] / lambdaT[i, i];
        for j in 1:nF-nR loop
          if i <> j then
           lambdaT[j,:] :=lambdaT[j, :] - lambdaT[j, i]*lambdaT[i, :]/lambdaT[i, i];
          end if;
        end for;
    else
      break;
    end if;
  end for;

  // get rid of negativ elements in nR+1:nF columns
  for i in nR+1:nF loop
    if min(lambdaT[:,i]) < 0 then
      for j in 1:nF-nR loop
        if lambdaT[j,i] < 0 then
          index :=Modelica.Math.BooleanVectors.firstTrueIndex({lambdaT[k,i] > 0 for k in 1:nF-nR});
          if index > 0 then
            lambdaT[j,:] :=lambdaT[j, :] + lambdaT[index, :] * abs(lambdaT[j,i])/lambdaT[index,i];
            index :=0;
          end if;
        end if;
      end for;
    end if;
  end for;

  // resort columns accordingly
  for i in 2:nF loop
    if resortindex[nF+1-i] then
      lambdaT[:,nF+1-i:nF] :=lambdaT[:, resort1[nF + 1 - i:nF,nF+1-i]];
    end if;
  end for;

  // transpose of lambdaT
  lambda :=transpose(lambdaT);

  //eliminate rounding errors
  for i in 1:size(lambda,1) loop
    for j in 1:size(lambda,2) loop
      if lambda[i,j] < 1e-15 then
        lambda[i,j] :=0;
      end if;
    end for;
  end for;

  else
    lambda :=identity(nF);
  end if;

  // deprecated
//    lambdaT :=total[nR + 1:nF, nR + 1:nR + nF];
//    for i in 1:nF-nR loop
//      if min(lambdaT) < 0 then
//        if not lambdaT[i,i] <> 0 then
//          (iDummy[i:nF-nR],sort[i:nF-nR]):=Modelica.Math.Vectors.sort(abs(lambdaT[i:nF-nR, i]),false);
//          for j in i:nF-nR loop
//            sort[j] :=sort[j] + (i-1);
//          end for;
//          lambdaT[i:nF-nR,:] :=lambdaT[sort[i:nF-nR], :];
//        end if;
//        lambdaT[i,:] :=lambdaT[i, :] / lambdaT[i, i];
//        for j in 1:nF-nR loop
//          if i <> j then
//            lambdaT[j,:] :=lambdaT[j, :] - lambdaT[j, i]*lambdaT[i, :]/lambdaT[i, i];
//          end if;
//        end for;
//      else
//        break;
//      end if;
//    end for;

//     for i in 1:nF-nR loop
//       if min(lambdaT) < 0 then
//          (iDummy[i:nF-nR],sort[i:nF-nR]):=Modelica.Math.Vectors.sort(lambdaT[i:nF-nR, i+nR],false);
//          for j in i:nF-nR loop
//            sort[j] :=sort[j] + (i-1);
//          end for;
//          lambdaT[i:nF-nR,:] :=lambdaT[sort[i:nF-nR], :];
//         for j in 1:nF-nR loop
//           if lambdaT[j,i+nR] < 0 then
//             lambdaT[j,:] :=lambdaT[j, :] - lambdaT[j, i+nR]*lambdaT[i, :]/lambdaT[i, i+nR];
//           end if;
//         end for;
//       else
//         break;
//       end if;
//     end for;
end calc_lambda;
