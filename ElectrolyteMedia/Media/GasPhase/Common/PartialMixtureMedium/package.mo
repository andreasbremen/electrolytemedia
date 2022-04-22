within ElectrolyteMedia.Media.GasPhase.Common;
partial package PartialMixtureMedium "Base class for pure substances of several chemical substances"
  extends Modelica.Media.Interfaces.PartialMedium(redeclare replaceable record
    FluidConstants =
        Modelica.Media.Interfaces.Types.IdealGas.FluidConstants);

  redeclare replaceable record extends ThermodynamicState
    "Thermodynamic state variables"
    Density d "density of medium";
    Temperature T "Temperature of medium";
    MassFraction[nX] X(start=reference_X)
      "Mass fractions (= (component mass)/total mass  m_i/m)";
  end ThermodynamicState;
  constant FluidConstants[nS] fluidConstants "Constant data for the fluid";

  replaceable function gasConstant
    "Return the gas constant of the mixture (also for liquids)"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "Thermodynamic state";
    output SI.SpecificHeatCapacity R "Mixture gas constant";
  end gasConstant;
end PartialMixtureMedium;
