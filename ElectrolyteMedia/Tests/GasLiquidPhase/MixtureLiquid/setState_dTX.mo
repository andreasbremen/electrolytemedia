within ElectrolyteMedia.Tests.GasLiquidPhase.MixtureLiquid;
model setState_dTX
  extends Base(
    redeclare package Medium =
        Media.GasLiquidPhase.MixtureLiquids.ExampleMedium,
    m0={0,0,0,1,1,0,0},
    Xg0={0,0.8,0.2},
    liquidFraction0=0.5);

  parameter SI.Density d0 = 3.1713522;

  Medium.ThermodynamicState state;
  Medium.BaseProperties medium(Xredstart = Xred0,Tstart=T0,pstart=p0);

equation

  state = Medium.setState_dTX(d0, T0, Xred0);

  medium.T = state.T;
  medium.p = state.p;
  medium.X[1:Medium.nX-1] = state.X[1:Medium.nX-1];

  annotation (experiment(Tolerance=1e-06));
end setState_dTX;
