within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_mfromY
  "Calculates molality vector from mole fraction vector of liquid species"
  input SI.MoleFraction[:] Y;
  output SI.MoleFraction[size(Y,1)-1] m;

algorithm
  for i in 1:size(Y,1)-1 loop
    m[i] :=Y[i]/IF97.MH2O/Y[end];//(max(1e-30,Y[i]))/IF97.MH2O/Y[end];//max(1e-30,Y[i]/IF97.MH2O/Y[end]);//max(Modelica.Constants.eps,Y[i]/IF97.MH2O/Y[end]);
  end for;
  annotation(Inline=true,smoothOrder=5);
end calc_mfromY;
