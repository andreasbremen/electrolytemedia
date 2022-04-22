within ElectrolyteMedia.Media.GasPhase.Common.PartialMixtureMedium;
function moleToMassFractions "Return mass fractions X from mole fractions"
  extends Modelica.Icons.Function;
  input SI.MoleFraction moleFractions[nX] "Mole fractions of mixture";
  input MolarMass[nX] MMX "Molar masses of components";
  output SI.MassFraction X[nX]
    "Mass fractions of gas mixture";
protected
  MolarMass Mmix=moleFractions*MMX "Molar mass of mixture";
algorithm
  for i in 1:nX loop
    X[i] := moleFractions[i]*MMX[i]/Mmix;
  end for;
  annotation (smoothOrder=5);
end moleToMassFractions;
