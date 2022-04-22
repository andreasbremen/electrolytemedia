within ElectrolyteMedia.Media.LiquidPhase.Common;
record DataRecord
  "Coefficient data record for properties of aqueous solute species"
  extends Modelica.Icons.Record;
//   parameter String name;
  parameter SI.MolarEnergy G_ref = 0 "Molar Gibbs free energy at standard state in J/mol";
  parameter SI.MolarEnergy H_ref = 0 "Molar enthalpy at reference state in J/mol";
  parameter SI.MolarEntropy S_ref = 0 "Molar entropy at reference state in J/mol/K";
  parameter Real a1 = 0 "HKF parameter in J/mol/Pa";
  parameter Real a2 = 0 "HKF parameter in J/mol";
  parameter Real a3 = 0 "HKF parameter in J/mol/K/Pa";
  parameter Real a4 = 0 "HKF parameter in J/mol/K";
  parameter Real c1 = 0 "HKF parameter in J/mol/K";
  parameter Real c2 = 0 "HKF parameter in J*K/mol";
  parameter Real w_ref = 0 "Conventional Born coefficient at standard state in J/mol";
  parameter Real z = 0 "charge of aqueous species";
  parameter Real a_0 = 0 "Kielland parameter for ion size in Debye Hückel equation in Angstrom";
  parameter SI.MolarMass MM = 0 "Molar mass of species";
  parameter SI.SpecificHeatCapacity R = 0 "Specific gas constant";
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end DataRecord;
