within ElectrolyteMedia.Tests.Solid;
partial model Base

//Declare Medium
replaceable package Medium =Media.SolidPhase.Common.MixtureSolid;

//Declare ThermodynamicState
parameter AbsolutePressure p_in=100e5;
parameter Temperature T_in=400;
Medium.ThermodynamicState state;

//Outputs that are calculated
AbsolutePressure p_calc;
Temperature T_calc;
 Density d_calc;
 SpecificEnthalpy h_calc;
 SpecificInternalEnergy u_calc;
 SpecificEntropy s_calc;
 SpecificEnergy g_calc;
   SpecificEnergy f_calc;

equation

p_calc = Medium.pressure(state);
T_calc = Medium.temperature(state);
 d_calc = Medium.density(state);
 h_calc = Medium.specificEnthalpy(state);
 u_calc = Medium.specificInternalEnergy(state);
 s_calc = Medium.specificEntropy(state);
 g_calc = Medium.specificGibbsEnergy(state);
 f_calc = Medium.specificHelmholtzEnergy(state);

end Base;
