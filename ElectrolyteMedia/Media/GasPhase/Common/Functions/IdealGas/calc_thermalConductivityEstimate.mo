within ElectrolyteMedia.Media.GasPhase.Common.Functions.IdealGas;
function calc_thermalConductivityEstimate
  "Thermal conductivity of polyatomic gases(Eucken and Modified Eucken correlation)"
  extends Modelica.Icons.Function;
  input Modelica.Media.Interfaces.Types.SpecificHeatCapacity Cp
    "Constant pressure heat capacity";
  input Modelica.Media.Interfaces.Types.DynamicViscosity eta
    "Dynamic viscosity";
  input Integer method(min=1,max=2)=1
    "1: Eucken Method, 2: Modified Eucken Method";
  output Modelica.Media.Interfaces.Types.ThermalConductivity lambda
    "Thermal conductivity [W/(m.k)]";
algorithm
  lambda := if method == 1 then eta*(Cp - datafun.R + (9/4)*datafun.R) else eta*(Cp
     - datafun.R)*(1.32 + 1.77/((Cp/Modelica.Constants.R) - 1.0));
  annotation (smoothOrder=2,
              Documentation(info="<html>
<p>
This function provides two similar methods for estimatinGfun the
thermal conductivity of polyatomic gases.
The Eucken method (input method == 1) gives good results for low temperatures,
but it tends to give an underestimated value of the thermal conductivity
(lambda) at higher temperatures.<br>
The Modified Eucken method (input method == 2) gives good results for
high-temperatures, but it tends to give an overestimated value of the
thermal conductivity (lambda) at low temperatures.
</p>
</html>"));
end calc_thermalConductivityEstimate;
