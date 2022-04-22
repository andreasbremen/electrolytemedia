within ElectrolyteMedia.Tests.Solid;
model setState_psX

  extends Base(redeclare package Medium =
        Media.SolidPhase.MixtureSolid.ExampleMedium);

  parameter SpecificEntropy s_in = 1371.9875;

equation
  state = Medium.setState_psX(p_in,s_in);

  annotation (experiment(Tolerance=1e-07));
end setState_psX;
