within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_alpha2 "Calculates alpha2 for salts of different valence types"
  output Real[nLifun,nLifun] alpha2;
algorithm
  for i in 1:nLifun loop
      for j in 1:nLifun loop
        if abs(datafun[i].z) >= 2 and abs(datafun[j].z) >= 2 then
          alpha2[i,j] :=12;
        else
          alpha2[i,j] :=0;
        end if;
      end for;
  end for;
end calc_alpha2;
