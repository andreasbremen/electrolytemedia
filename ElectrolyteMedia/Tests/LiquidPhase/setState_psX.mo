within ElectrolyteMedia.Tests.LiquidPhase;
model setState_psX
  extends Base(
    redeclare package Medium =
        Media.LiquidPhase.MixtureLiquids.ExampleMedium,
    m_in = {5,0,4,0,4,0});

  parameter SI.SpecificEntropy s_in = 4704.31;//3400.0916;

  Medium.BaseProperties medium(Xredstart=Xred_in);
equation
  medium.T=state.T;
  medium.p=state.p;
  medium.Xi=state.X[1:Medium.nX-1];
  state = Medium.setState_psX(p_in, s_in,Xred_in);

  annotation (experiment(Tolerance=1e-06));
end setState_psX;
