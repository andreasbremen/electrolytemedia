within ElectrolyteMedia.Tests.LiquidPhase;
model setState_phX
  extends LiquidPhase.Base(
    redeclare package Medium =
        Media.LiquidPhase.MixtureLiquids.ExampleMedium,
    m_in = {5,0,4,0,4,0});

  parameter SI.SpecificEnthalpy h_in = -13957863;

  Medium.BaseProperties medium(Xredstart=Xred_in);
equation
  medium.T=state.T;
  medium.p=state.p;
  medium.Xi=state.X[1:Medium.nX-1];
  state = Medium.setState_phX(p_in, h_in,Xred_in);

  annotation (experiment(Tolerance=1e-06));
end setState_phX;
