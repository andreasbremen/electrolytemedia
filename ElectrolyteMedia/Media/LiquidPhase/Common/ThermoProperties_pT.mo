within ElectrolyteMedia.Media.LiquidPhase.Common;
record ThermoProperties_pT "Thermodynamic property data for pressure p and temperature T as dynamic states"
  extends Modelica.Media.Common.ThermoFluidSpecial.ThermoProperties_pT(d(min=0));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ThermoProperties_pT;
