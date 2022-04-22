within ElectrolyteMedia.Tests.GasLiquidPhase.MixtureLiquid;
model Test_volume_Modelica
  extends Modelica.Icons.Example;
  extends Modelica.Media.Examples.Utilities.PartialTestModel(
      redeclare package Medium =
        Media.GasLiquidPhase.MixtureLiquids.ExampleMedium,
    fixedMassFlowRate(
      use_T_ambient=true,
      T_ambient=323.15,
      X_ambient={0.32131994,0.0090731615,0.3474194,0.3216091,5.78369e-15,0.000578369}),
    volume(V=1));
  annotation (experiment(
      StopTime=50,
      Tolerance=1e-08,
      __Dymola_Algorithm="Lsodar"));
end Test_volume_Modelica;
