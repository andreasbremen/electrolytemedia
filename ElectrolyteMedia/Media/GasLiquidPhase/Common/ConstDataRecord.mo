within ElectrolyteMedia.Media.GasLiquidPhase.Common;
record ConstDataRecord
  constant Real eta = 1.66027e5 "constant for w calculations";
  constant Modelica.SIunits.Pressure Psi = 2600e5 "solvent characerstic constant in Pa";
  constant Modelica.SIunits.Temperature Theta = 228 "solvent characteristic constant in K";
  constant Real Z_ref = -1.278034682000000E-002 "Z at reference state";
  constant Real Y_ref = -5.798650444000000E-005 "Y at reference state";
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ConstDataRecord;
