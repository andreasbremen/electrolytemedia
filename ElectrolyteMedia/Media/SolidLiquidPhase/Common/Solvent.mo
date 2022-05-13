within ElectrolyteMedia.Media.SolidLiquidPhase.Common;
record Solvent
  constant Real a[:] = {-2.037662, 5.747e-3, -6.557892e-6};
  constant Real b[:] = {6.107361, -1.074377e-2, 1.268348e-5};
  constant Modelica.SIunits.Density rho_0 = 1000;
  constant Real a_f[:] = {3.66666e-16, -1.504956e-10, 5.01799e-14};
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Solvent;
