within ElectrolyteMedia.Tests.Solid;
model setState_phX

  extends Base(redeclare package Medium =
        Media.SolidPhase.MixtureSolid.ExampleMedium);

  parameter SpecificEnthalpy h_in = -6987912;//-14539058;

equation
  state = Medium.setState_phX(p_in, h_in);

  annotation (experiment(Tolerance=1e-07));
end setState_phX;
