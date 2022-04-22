within ElectrolyteMedia.Tests.GasLiquidPhase.SaturatedLiquid;
model setState_phX
  extends SaturatedLiquid.Base(
    redeclare package Medium =
        Media.GasLiquidPhase.SaturatedLiquids.ExampleMedium,
    m0={0,0,0.2523,0},
    Xg0={1},
    liquidFraction0=1,
    p0=1e5);  //0.969702e5);

  parameter SI.SpecificEnthalpy h0 = -15491761;

  Medium.ThermodynamicState state;
  Medium.BaseProperties medium(Xredstart = Xred0,Tstart=T0,pstart=p0);

  SpecificEntropy s_calc;
  SpecificEnthalpy h_calc;
  Density d_calc;
equation

  s_calc = Medium.specificEntropy(state);
  h_calc = Medium.specificEnthalpy(state);
  d_calc = Medium.density(state);

  state = Medium.setState_phX(p0, h0, Xred0);

  medium.T = state.T;
  medium.p = state.p;
  medium.X[1:Medium.nX-1] = state.X[1:Medium.nX-1];

  annotation (experiment(Tolerance=1e-06));
end setState_phX;
