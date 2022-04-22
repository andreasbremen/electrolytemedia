within ElectrolyteMedia.Media.LiquidPhase;
package Common "Contains functions and models to calculate ideal gas thermodynamics"

  constant Modelica.SIunits.Pressure pref = 1e5 "Pressure at reference point";
  constant Modelica.SIunits.Temperature Tref = 298.15 "Temperature at reference point";
  constant Real eta = 1.66027e5 "Constant for w calculations";
  constant Modelica.SIunits.Pressure Psi = 2600e5 "Solvent characerstic constant in Pa";
  constant Modelica.SIunits.Temperature Theta = 228 "Solvent characteristic constant in K";
  constant Real Zref = -1.278034682000000E-002 "Z at reference state";
  constant Real Yref = -5.798650444000000E-005 "Y at reference state";



end Common;
