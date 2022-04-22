within ElectrolyteMedia.Tests.LiquidPhase;
model ModelicaExample
  extends Modelica.Media.Examples.Utilities.PartialTestModel(
         redeclare package Medium =
        Media.LiquidPhase.MixtureLiquids.ExampleMedium, volume(
      T_start=298.15),
    fixedMassFlowRate(X_ambient=X_start));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(Tolerance=1e-08, __Dymola_Algorithm="Dassl"));
end ModelicaExample;
