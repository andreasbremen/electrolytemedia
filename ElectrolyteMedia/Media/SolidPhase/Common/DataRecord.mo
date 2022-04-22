within ElectrolyteMedia.Media.SolidPhase.Common;
record DataRecord "Coefficient data record for properties of solid species"

//   parameter Real[SolidModel2.nX] X_default=fill(1/SolidModel2.nX,
//       SolidModel2.nX);
    parameter Modelica.SIunits.SpecificHeatCapacity R;
//   parameter Modelica.SIunits.Temperature T_0 = 298.15 "Reference temperature";
//   parameter Modelica.SIunits.Pressure p_0 = 1e5 "Reference pressure";
  parameter Modelica.SIunits.MolarEnergy G_ref = 0 "Molar Gibbs free energy at standard state in J/mol";
  parameter Modelica.SIunits.MolarEnthalpy H_ref = 0 "Molar enthalpy at reference state in J/mol";
  parameter Modelica.SIunits.MolarEntropy S_ref = 0 "Molar entropy at reference state in J/mol/K";
  parameter Modelica.SIunits.MolarVolume V_ref = 0 "Molar volume at reference state in m3/mol";
  parameter Real a = 0 "coefficient of c_p polynomial";
  parameter Real b = 0 "coefficient of c_p polynomial";
  parameter Real c = 0 "coefficient of c_p polynomial";
  parameter Real d = 0 "coefficient of c_p polynomial";
  parameter Real alpha_0 = 0 "parameter for solid EOS";
  parameter Real k_0 = 0 "parameter for solid EOS";
  parameter Real k_0_ = 0 "parameter for solid EOS";
  parameter Real k_0__ = 0 "parameter for solid EOS";
  parameter Real n = 0 "number of atoms";
  parameter Modelica.SIunits.MolarMass MM "Molar mass of mineral";

//    parameter String name;
//   parameter Modelica.SIunits.MolarMass MM;

//    parameter String name1;
//    parameter String name2;
//    parameter String name3;
//    parameter String name4;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end DataRecord;
