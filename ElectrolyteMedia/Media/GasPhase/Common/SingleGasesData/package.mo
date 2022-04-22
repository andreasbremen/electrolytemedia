within ElectrolyteMedia.Media.GasPhase.Common;
package SingleGasesData "Properties of fluids from SUPCRTBL by Zimmer(2016) and Johnson(1992)"

  constant GasDataRecord[:] H2O_CO2_O2 = {H2O,CO2,O2} "Water,carbon dioxide,oxygen";
  constant GasDataRecord[:] H2O_CO2 = {H2O,CO2} "Water,carbon dioxide";
  constant GasDataRecord[:] CO2_H2O = {CO2,H2O} "Carbon dioxide,water";
  constant GasDataRecord[:] H2O_O2 = {H2O,O2} "Water,oxygen";
  constant GasDataRecord[:] H2O_array = {H2O} "Water";
protected 
  constant GasDataRecord CH4(
    MM=1.604160e-02,
    H_ref=-7.481000e+04,
    G_ref=-5.071000e+04,
    S_ref=1.862600e+02,
    a=1.501000e+02,
    b=2.063000e-03,
    c=3.427700e+06,
    d=-2.650400e+03,
    R = Modelica.Constants.R/1.604160e-2,
    T_c=FluidData.CH4.criticalTemperature,
    p_c=FluidData.CH4.criticalPressure,
    w=FluidData.CH4.acentricFactor) "Methane";
//     name="CH4",

  constant GasDataRecord CO(
    MM=2.801000e-02,
    H_ref=-1.105300e+05,
    G_ref=-1.371700e+05,
    S_ref=1.976700e+02,
    a=4.570000e+01,
    b=-9.700000e-05,
    c=6.627000e+05,
    d=-4.147000e+02,
    R = Modelica.Constants.R/2.801000e-02,
    T_c=FluidData.CO.criticalTemperature,
    p_c=FluidData.CO.criticalPressure,
    w=FluidData.CO.acentricFactor) "Carbon monoxide";
//     name="CO",

  constant GasDataRecord CO2(
    MM=4.401000e-02,
    H_ref=-3.935100e+05,
    G_ref=-3.943500e+05,
    S_ref=2.137000e+02,
    a=8.780000e+01,
    b=-2.644000e-03,
    c=7.064000e+05,
    d=-9.989000e+02,
    R = Modelica.Constants.R/4.401000e-02,
    T_c=FluidData.CO2.criticalTemperature,
    p_c=FluidData.CO2.criticalPressure,
    w=FluidData.CO2.acentricFactor) "Carbon dioxide";
//     name="CO2",

  constant GasDataRecord H2(
    MM=2.015800e-03,
    H_ref=0.000000e+00,
    G_ref=-1.000000e+01,
    S_ref=1.307000e+02,
    a=2.330000e+01,
    b=4.627000e-03,
    c=0.000000e+00,
    d=7.630000e+01,
    R = Modelica.Constants.R/2.015800e-03,
    T_c=FluidData.H2.criticalTemperature,
    p_c=FluidData.H2.criticalPressure,
    w=FluidData.H2.acentricFactor) "Hydrogen";
//     name="H2",

  constant GasDataRecord H2S(
    MM=3.407580e-02,
    H_ref=-2.030000e+04,
    G_ref=-3.313000e+04,
    S_ref=2.057700e+02,
    a=4.740000e+01,
    b=1.024000e-02,
    c=6.159000e+05,
    d=-3.978000e+02,
    R = Modelica.Constants.R/3.407580e-02,
    T_c=FluidData.H2S.criticalTemperature,
    p_c=FluidData.H2S.criticalPressure,
    w=FluidData.H2S.acentricFactor) "Hydrogen sulfide";
//     name="H2S",

  constant GasDataRecord O2(
    MM=3.200000e-02,
    H_ref=0.000000e+00,
    G_ref=-1.000000e+01,
    S_ref=2.052000e+02,
    a=4.830000e+01,
    b=-6.910000e-04,
    c=4.992000e+05,
    d=-4.207000e+02,
    R = Modelica.Constants.R/3.200000e-02,
    T_c=FluidData.O2.criticalTemperature,
    p_c=FluidData.O2.criticalPressure,
    w=FluidData.O2.acentricFactor) "Oxygen";
//     name="O2",

  constant GasDataRecord S2(
    MM=6.412000e-02,
    H_ref=1.285400e+05,
    G_ref=7.878000e+04,
    S_ref=2.310000e+02,
    a=3.710000e+01,
    b=2.398000e-03,
    c=-1.610000e+05,
    d=-6.500000e+01,
    R = Modelica.Constants.R/6.412000e-02,
    T_c=FluidData.S2.criticalTemperature,
    p_c=FluidData.S2.criticalPressure,
    w=FluidData.S2.acentricFactor) "Sulfur";
//     name="S2",

  constant GasDataRecord H2O(
    MM=1.80160e-02,
    H_ref=-2.418100e+05,
    G_ref=-2.285600e+05,
    S_ref=1.888000e+02,
    a=4.010000e+01,
    b=8.656000e-03,
    c=4.875000e+05,
    d=-2.512000e+02,
    R = Modelica.Constants.R/1.80160e-02,
    T_c=FluidData.H2O.criticalTemperature,
    p_c=FluidData.H2O.criticalPressure,
    w=FluidData.H2O.acentricFactor) "Water";
    //     name="H2O",

end SingleGasesData;
