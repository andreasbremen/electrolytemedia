within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_gamma_ideal
  "Calculates activity coefficient of solutes and solvent at T and p for ideal case"
  output SI.SpecificEnergy[nLfun] gamma;
algorithm
    for i in 1:nLfun loop
      if i<nLfun then
        gamma[i] :=1;
      else
        gamma[i] :=MixtureLiquid.MH2O;
      end if;
    end for;

  annotation(smoothOrder=5);
end calc_gamma_ideal;
