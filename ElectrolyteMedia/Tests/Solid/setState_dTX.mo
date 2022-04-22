within ElectrolyteMedia.Tests.Solid;
model setState_dTX

  extends Base(redeclare package Medium =
        Media.SolidPhase.MixtureSolid.ExampleMedium);

  parameter Density d_in = 2081.3136;

equation
  state = Medium.setState_dTX(d_in, 350);

  annotation (experiment(Tolerance=1e-07));
end setState_dTX;
