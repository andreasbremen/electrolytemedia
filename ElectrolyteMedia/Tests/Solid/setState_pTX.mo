within ElectrolyteMedia.Tests.Solid;
model setState_pTX

  extends Base(redeclare package Medium =
        Media.SolidPhase.MixtureSolid.ExampleMedium);

  MolarMass MM;

equation

  state = Medium.setState_pTX(p_in, 350);
  MM = Medium.molarMass(state);

end setState_pTX;
