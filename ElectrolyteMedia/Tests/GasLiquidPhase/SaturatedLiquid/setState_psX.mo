within ElectrolyteMedia.Tests.GasLiquidPhase.SaturatedLiquid;
model setState_psX
  extends SaturatedLiquid.Base(
    redeclare package Medium =
        Media.GasLiquidPhase.SaturatedLiquids.ExampleMedium,
    m0={0,0,0.2532,0},
    Xg0={1},
    liquidFraction0=1,
    p0=1e5);

  parameter SI.SpecificEntropy s0 = 4842.8145;

  Medium.ThermodynamicState state;
  Medium.BaseProperties medium(Xredstart = Xred0,pstart=p0,Tstart=T0);

  SpecificEntropy s_calc;
  SpecificEnthalpy h_calc;
  Density d_calc;
equation

  s_calc = Medium.specificEntropy(state);
  h_calc = Medium.specificEnthalpy(state);
  d_calc = Medium.density(state);

  state = Medium.setState_psX(p0, s0, Xred0);

  medium.T = state.T;
  medium.p = state.p;
  medium.X[1:Medium.nX-1] = state.X[1:Medium.nX-1];

  annotation (experiment(Tolerance=1e-06));
end setState_psX;
