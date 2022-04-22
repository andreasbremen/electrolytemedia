within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_alpha1 "Calculates alpha for salts of different valence types"
  output Real[nLifun,nLifun] alpha1;
algorithm
  for i in 1:nLifun loop
      for j in 1:nLifun loop
        if abs(datafun[i].z) >= 2 and abs(datafun[j].z) >= 2 then
          alpha1[i,j] :=1.4;
        elseif abs(datafun[i].z) == 0 or abs(datafun[j].z) == 0 then
          alpha1[i,j] :=0;
        else
          alpha1[i,j] :=2;
        end if;
      end for;
  end for;
end calc_alpha1;
