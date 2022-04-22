within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid;
function molarVolume "Calculates molar volume of mixture"
  input ThermodynamicState state;
  output MolarVolume v;
protected
  Density d=Functions.calc_d(
        state.T,
        state.p,
        state.X);
algorithm
  v := molarMass(state)/d;
  annotation(Inline = true, smoothOrder = 3);
end molarVolume;
