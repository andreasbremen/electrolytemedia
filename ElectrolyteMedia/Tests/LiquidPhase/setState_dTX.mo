within ElectrolyteMedia.Tests.LiquidPhase;
model setState_dTX
  extends Base(
    redeclare package Medium =
        Media.LiquidPhase.MixtureLiquids.ExampleMedium,
    m_in = {5,0,4,0,4,0});

  parameter SI.Density d_in = 1.1312762e3;

  Medium.BaseProperties medium(Xredstart=Xred_in);
equation
  medium.T=state.T;
  medium.p=state.p;
  medium.Xi=state.X[1:Medium.nX-1];
  state = Medium.setState_dTX(d_in,T_in,Xred_in);

  annotation (experiment(Tolerance=1e-06));
end setState_dTX;
