within ElectrolyteMedia.Tests.GasLiquidPhase.MixtureLiquid;
partial model Base "Test method to generate outputs from given inputs"

  //Declare Medium
  replaceable package Medium =
      Media.GasLiquidPhase.Common.MixtureLiquid;

  //Declare ThermodynamicState
  //To be specified
  parameter Real[Medium.nL-1] m0;
  parameter SI.MassFraction[Medium.nG] Xg0;
  parameter Real liquidFraction0;
  //Given
  parameter Temperature T0=300;
  parameter AbsolutePressure p0=1e5;
  parameter SI.MassFraction[Medium.nL] Xl0 = Medium.Functions.LiquidFunctions.calc_Xfromm(m0);
  parameter SI.MassFraction[Medium.nF] X0 = cat(1,(1-liquidFraction0)*Xg0,liquidFraction0*Xl0);
  parameter SI.MassFraction[Medium.nX] Xred0 = transpose(Medium.lambda_mass)*X0/sum(transpose(Medium.lambda_mass)*X0);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(Tolerance=1e-08));
end Base;
