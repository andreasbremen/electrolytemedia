within ElectrolyteMedia.Tests.GasPhase;
model ModelicaTestModel_PengRobinson
  extends Modelica.Media.Examples.Utilities.PartialTestModel(
         redeclare package Medium =
        Media.GasPhase.MixtureGases.CO2_H2O_PR);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=10,
      Tolerance=1e-06,
      __Dymola_Algorithm="Dassl"));
end ModelicaTestModel_PengRobinson;
