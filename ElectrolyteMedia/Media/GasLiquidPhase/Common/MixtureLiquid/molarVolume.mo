within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid;
function molarVolume "Calculates molar volume of mixture"
  input ThermodynamicState state;
  output MolarVolume v;
protected
  Density d = density(state);
algorithm
  v := molarMass(state)/d;
  annotation(Inline = true, smoothOrder = 3);
end molarVolume;
