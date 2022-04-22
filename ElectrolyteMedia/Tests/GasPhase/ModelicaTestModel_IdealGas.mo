within ElectrolyteMedia.Tests.GasPhase;
model ModelicaTestModel_IdealGas
  extends Modelica.Media.Examples.Utilities.PartialTestModel(redeclare package
      Medium = Media.GasPhase.MixtureGases.CO2_H2O_IG);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end ModelicaTestModel_IdealGas;
