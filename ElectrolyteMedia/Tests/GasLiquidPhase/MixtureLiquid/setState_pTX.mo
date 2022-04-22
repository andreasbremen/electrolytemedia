within ElectrolyteMedia.Tests.GasLiquidPhase.MixtureLiquid;
model setState_pTX
  extends Base(
    redeclare package Medium =
        Media.GasLiquidPhase.MixtureLiquids.ExampleMedium,
    m0={0,0,0,1,1,0,0},
    Xg0={0,0.8,0.2},
    liquidFraction0=0.5);

  Medium.ThermodynamicState state;
  Medium.BaseProperties medium(Xredstart = Xred0,Tstart=T0,pstart=p0);

  SpecificEntropy s_calc;
  SpecificEnthalpy h_calc;
  Density d_calc;

equation

  s_calc = Medium.specificEntropy(state);
  h_calc = Medium.specificEnthalpy(state);
  d_calc = Medium.density(state);

  state = Medium.setState_pTX(p0, T0, Xred0);

  medium.T = state.T;
  medium.p = state.p;
  medium.X[1:Medium.nX-1] = state.X[1:Medium.nX-1];

  annotation (experiment(Tolerance=1e-06));
end setState_pTX;
