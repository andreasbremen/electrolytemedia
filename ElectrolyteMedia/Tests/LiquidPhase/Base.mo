within ElectrolyteMedia.Tests.LiquidPhase;
partial model Base
  "Test method to generate outputs from given inputs"

  //Declare Medium
  replaceable package Medium =Media.LiquidPhase.Common.MixtureLiquid;

  //Declare ThermodynamicState
  parameter Temperature T_in=300;
  parameter AbsolutePressure p_in=10e5;
  parameter Real[Medium.nF-1] m_in;
  parameter SI.MassFraction[Medium.nF] X_in = Medium.Functions.calc_Xfromm(m_in);
  parameter SI.MassFraction[Medium.nX] Xred_in = transpose(Medium.lambda_mass)*X_in/sum(transpose(Medium.lambda_mass)*X_in);
  Medium.ThermodynamicState state;

  //Outputs that are calculated
  Temperature T_calc;
  Density d_calc;
  AbsolutePressure p_calc;
  SpecificHeatCapacity cp_calc;
  SpecificEnthalpy h_calc;
  SpecificInternalEnergy u_calc;
  SpecificEntropy s_calc;
  SpecificEntropy R_calc;
  MoleFraction[Medium.nF] Yfull;

equation

  T_calc = Medium.temperature(state);
  d_calc = Medium.density(state);
  p_calc = Medium.pressure(state);
  cp_calc = Medium.heatCapacity_cp(state);
  h_calc = Medium.specificEnthalpy(state);
  u_calc = Medium.specificInternalEnergy(state);
  s_calc = Medium.specificEntropy(state);
  R_calc = Medium.gasConstant(state);

  Yfull = Medium.massToMoleFractions(state.Xfull,Medium.MMX);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(Tolerance=1e-08));
end Base;
