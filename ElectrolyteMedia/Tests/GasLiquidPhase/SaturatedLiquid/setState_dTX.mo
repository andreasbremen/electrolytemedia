within ElectrolyteMedia.Tests.GasLiquidPhase.SaturatedLiquid;
model setState_dTX
  extends SaturatedLiquid.Base(
    redeclare package Medium =
        Media.GasLiquidPhase.SaturatedLiquids.ExampleMedium,
    m0={0,0,0.2352,0},
    Xg0={1},
    liquidFraction0=1,
    T0=101.10931 + 273.15);

  parameter SI.Density d0 = 67.57274;

  Medium.ThermodynamicState state;
  Medium.BaseProperties medium(Xredstart = Xred0,Tstart=T0,pstart=p0);

equation

  state = Medium.setState_dTX(d0, T0, Xred0);

  medium.T = state.T;
  medium.p = state.p;
  medium.X[1:Medium.nX-1] = state.X[1:Medium.nX-1];

  annotation (experiment(Tolerance=1e-06));
end setState_dTX;
