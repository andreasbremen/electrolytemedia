within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_Yfromm
  "Calculates mole fraction vector from molality vector of liquid species"
  input Real[nLifun] m;
  output SI.MoleFraction[nLifun+1] Y;

algorithm
  Y[nLifun+1] :=1/IF97.MH2O/(sum(m) + 1/IF97.MH2O);
  for i in 1:nLifun loop
    Y[i] :=m[i]*IF97.MH2O*Y[nLifun+1];
  end for;

  annotation(smoothOrder=5);
end calc_Yfromm;
