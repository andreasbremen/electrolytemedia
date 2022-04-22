within ElectrolyteMedia.Media.GasPhase.Common;
package FluidData "Properties of fluids from SUPCRTBL by Zimmer(2016) and Johnson(1992)"

  constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants[:] CO2_H2O = {CO2,H2O} "Carbon dioxide,water";

  constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants CH4(
                       chemicalFormula =        "CH4",
                       iupacName =              "unknown",
                       structureFormula =       "unknown",
                       casRegistryNumber =      "74-82-8",
                       meltingPoint =            90.69,
                       normalBoilingPoint =     111.66,
                       criticalTemperature =    190.56,
                       criticalPressure =        45.99e5,
                       criticalMolarVolume =     98.60e-6,
                       acentricFactor =           0.01,
                       dipoleMoment =             0.0,
                       molarMass =              1.604160e-02,
                       hasDipoleMoment =       true,
                       hasIdealGasHeatCapacity=true,
                       hasCriticalData =       true,
                       hasAcentricFactor =     true);

  constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants CO(
                       chemicalFormula =        "CO",
                       iupacName =              "unknown",
                       structureFormula =       "unknown",
                       casRegistryNumber =      "630-08-0",
                       meltingPoint =            68.15,
                       normalBoilingPoint =      81.66,
                       criticalTemperature =    132.92,
                       criticalPressure =        34.99e5,
                       criticalMolarVolume =     93.10e-6,
                       acentricFactor =           0.093,
                       dipoleMoment =             0.1,
                       molarMass =              2.801000e-02,
                       hasDipoleMoment =       true,
                       hasIdealGasHeatCapacity=true,
                       hasCriticalData =       true,
                       hasAcentricFactor =     true);

  constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants CO2(
                       chemicalFormula =        "CO2",
                       iupacName =              "unknown",
                       structureFormula =       "unknown",
                       casRegistryNumber =      "124-38-9",
                       meltingPoint =           216.58,
                       normalBoilingPoint =     -1.0,
                       criticalTemperature =    304.21,
                       criticalPressure =        73.83e5,
                       criticalMolarVolume =     94.07e-6,
                       acentricFactor =           0.231,
                       dipoleMoment =             0.0,
                       molarMass =              4.401000e-02,
                       hasDipoleMoment =       true,
                       hasIdealGasHeatCapacity=true,
                       hasCriticalData =       true,
                       hasAcentricFactor =     true);

  constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants H2(
                       chemicalFormula =        "H2",
                       iupacName =              "unknown",
                       structureFormula =       "unknown",
                       casRegistryNumber =      "800000-51-5",
                       meltingPoint =            13.56,
                       normalBoilingPoint =      20.38,
                       criticalTemperature =     33.19,
                       criticalPressure =        13.13e5,
                       criticalMolarVolume =     65.00e-6,
                       acentricFactor =          -0.216,
                       dipoleMoment =             0.0,
                       molarMass =              2.015800e-03,
                       hasDipoleMoment =       true,
                       hasIdealGasHeatCapacity=true,
                       hasCriticalData =       true,
                       hasAcentricFactor =     true);

  constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants H2S(
                       chemicalFormula =        "H2S",
                       iupacName =              "unknown",
                       structureFormula =       "unknown",
                       casRegistryNumber =      "7783-06-4",
                       meltingPoint =            187.45,
                       normalBoilingPoint =      212.95,
                       criticalTemperature =     373.53,
                       criticalPressure =        89.6291e5,
                       criticalMolarVolume =     98.00e-6,
                       acentricFactor =          0.1,
                       dipoleMoment =             0.0,
                       molarMass =              3.407580e-02,
                       hasDipoleMoment =       true,
                       hasIdealGasHeatCapacity=true,
                       hasCriticalData =       true,
                       hasAcentricFactor =     true);

  constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants O2(
                       chemicalFormula =        "O2",
                       iupacName =              "unknown",
                       structureFormula =       "unknown",
                       casRegistryNumber =      "7782-44-7",
                       meltingPoint =            54.36,
                       normalBoilingPoint =      90.17,
                       criticalTemperature =    154.58,
                       criticalPressure =        50.43e5,
                       criticalMolarVolume =     73.37e-6,
                       acentricFactor =         0.019,
                       dipoleMoment =           0.0,
                       molarMass =              3.200000e-02,
                       hasDipoleMoment =       true,
                       hasIdealGasHeatCapacity=true,
                       hasCriticalData =       true,
                       hasAcentricFactor =     true);

  constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants S2(
                       chemicalFormula =        "S2",
                       iupacName =              "unknown",
                       structureFormula =       "unknown",
                       casRegistryNumber =      "23550-45-0",
                       meltingPoint =            1,
                       normalBoilingPoint =      1,
                       criticalTemperature =    298.15,
                       criticalPressure =        1e5,
                       criticalMolarVolume =     1e-6,
                       acentricFactor =         0,
                       dipoleMoment =           0.0,
                       molarMass =              6.412000e-02,
                       hasDipoleMoment =       true,
                       hasIdealGasHeatCapacity=true,
                       hasCriticalData =       true,
                       hasAcentricFactor =     true);

  constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants H2O(
                       chemicalFormula =        "H2O",
                       iupacName =              "oxidane",
                       structureFormula =       "H2O",
                       casRegistryNumber =      "7732-18-5",
                       meltingPoint =           273.15,
                       normalBoilingPoint =     373.124,
                       criticalTemperature =    647.096,
                       criticalPressure =       220.64e5,
                       criticalMolarVolume =     55.95e-6,
                       acentricFactor =           0.348,
                       dipoleMoment =             1.8,
                       molarMass =              1.80160e-02,
                       hasDipoleMoment =       true,
                       hasIdealGasHeatCapacity=true,
                       hasCriticalData =       true,
                       hasAcentricFactor =     true);

end FluidData;
