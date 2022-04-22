within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp;
function calc_eta
  "Calculates dynamic viscosity of solvent water with IF97 formulation"
  input SI.Pressure p;
  input SI.Temperature T;
  output SI.DynamicViscosity eta;
protected
  SI.Density d;
algorithm
  d :=calc_rho(T, p);
  eta :=Modelica.Media.Water.IF97_Utilities.BaseIF97.Transport.visc_dTp(
    d,
    T,
    p);
end calc_eta;
