within ElectrolyteMedia.Media.GasPhase.Common;
record GasDataRecord
  "Data record containing coefficients for Peng-Robinson gasesCoefficient data record for properties of ideal gases"
  extends Modelica.Icons.Record;
  SI.MolarMass MM = 0 "Molar mass";
  SI.MolarEnthalpy H_ref = 0 "Enthalpy of formation at 298.15K";
  SI.MolarEnergy G_ref = 0 "Molar Gibbs free energy at standard state in J/mol";
  SI.MolarEntropy S_ref = 0 "Molar entropy at reference state in J/mol/K";
  Real a = 0 "Temperature coefficient";
  Real b = 0 "Temperature coefficient";
  Real c = 0 "Temperature coefficient";
  Real d = 0 "Temperature coefficient";
  SI.SpecificHeatCapacity R = 0 "Gas constant";
  SI.Temperature T_c = 0 "critical Temperature of species";
  SI.Pressure p_c = 0 "critical pressure of species";
  Real w = 0 "acentric factor of species";

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end GasDataRecord;
