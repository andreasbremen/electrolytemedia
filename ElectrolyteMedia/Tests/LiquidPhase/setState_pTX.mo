within ElectrolyteMedia.Tests.LiquidPhase;
model setState_pTX
  extends Base(
    redeclare package Medium =
        Media.LiquidPhase.MixtureLiquids.ExampleMedium,
    m_in = {5,0,4,0,4,0});

  Medium.BaseProperties medium(Xredstart=Xred_in);

equation
  medium.T=state.T;
  medium.p=state.p;
  medium.Xi=state.X[1:Medium.nX-1];
  state = Medium.setState_pTX(p_in, T_in, Xred_in);
  annotation (experiment(Tolerance=1e-08));
end setState_pTX;
