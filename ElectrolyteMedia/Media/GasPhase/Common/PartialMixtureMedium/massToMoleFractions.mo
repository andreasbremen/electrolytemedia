within ElectrolyteMedia.Media.GasPhase.Common.PartialMixtureMedium;
function massToMoleFractions "Return mole fractions from mass fractions X"
  extends Modelica.Icons.Function;
  input SI.MassFraction X[nX] "Mass fractions of mixture";
  input SI.MolarMass[nX] MMX "Molar masses of components";
  output SI.MoleFraction moleFractions[nX]
    "Mole fractions of gas mixture";
protected
  Real invMMX[nX] "Inverses of molar weights";
  SI.MolarMass Mmix "Molar mass of mixture";
algorithm
  for i in 1:nX loop
    invMMX[i] := 1/MMX[i];
  end for;
  Mmix := 1/(X*invMMX);
  for i in 1:nX loop
    moleFractions[i] := Mmix*X[i]/MMX[i];
  end for;
  annotation (smoothOrder=5);
end massToMoleFractions;
