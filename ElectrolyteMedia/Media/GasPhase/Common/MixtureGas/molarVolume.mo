within ElectrolyteMedia.Media.GasPhase.Common.MixtureGas;
function molarVolume "Calculates molar volume of mixture"
  input ThermodynamicState state;
  output MolarVolume v;
algorithm
  v := molarMass(state)/state.d;
  annotation(Inline = true, smoothOrder = 3);
end molarVolume;
