within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp.Reduced;
function calc_g "Gibbs function for region 1: g(p,T)"
  input SI.Temperature T "Temperature (K)";
  input SI.Pressure p "Pressure";
  output Real g
    "Dimensionless Gibbs function and derivatives w.r.t. pi and tau";
protected
  Real pi1 "Dimensionless pressure";
  Real tau1 "Dimensionless temperature";
  Real[45] o "Vector of auxiliary variables";
  Real pl "Auxiliary variable";
  Real pi(unit="1") "Dimensionless pressure";
  Real tau(unit="1") "Dimensionless temperature";
algorithm
  //         pl := min(p, IF97.PCRIT - 1);
//   assert(p > IF97.ptriple, "IF97 medium function g1 called with too low pressure\n"
//      + "p = " + String(p) + " Pa <= " + String(IF97.ptriple)
//      + " Pa (triple point pressure)");
//   assert(p <= 100.0e6, "IF97 medium function g1: the input pressure (= "
//      + String(p) + " Pa) is higher than 100 Mpa");
//   assert(T >= 273.15, "IF97 medium function g1: the temperature (= " +
//     String(T) + " K) is lower than 273.15 K!");
  pi := p/IF97.PSTAR1;
  tau := IF97.TSTAR1/T;
  pi1 := 7.1000000000000 - pi;
  tau1 := -1.22200000000000 + tau;
  o[1] := tau1*tau1;
  o[2] := o[1]*o[1];
  o[3] := o[2]*o[2];
  o[4] := o[3]*tau1;
  o[5] := 1/o[4];
  o[6] := o[1]*o[2];
  o[7] := o[1]*tau1;
  o[8] := 1/o[7];
  o[9] := o[1]*o[2]*o[3];
  o[10] := 1/o[2];
  o[11] := o[2]*tau1;
  o[12] := 1/o[11];
  o[13] := o[2]*o[3];
  o[14] := 1/o[3];
  o[15] := pi1*pi1;
  o[16] := o[15]*pi1;
  o[17] := o[15]*o[15];
  o[18] := o[17]*o[17];
  o[19] := o[17]*o[18]*pi1;
  o[20] := o[15]*o[17];
  o[21] := o[3]*o[3];
  o[22] := o[21]*o[21];
  o[23] := o[22]*o[3]*tau1;
  o[24] := 1/o[23];
  o[25] := o[22]*o[3];
  o[26] := 1/o[25];
  o[27] := o[1]*o[2]*o[22]*tau1;
  o[28] := 1/o[27];
  o[29] := o[1]*o[2]*o[22];
  o[30] := 1/o[29];
  o[31] := o[1]*o[2]*o[21]*o[3]*tau1;
  o[32] := 1/o[31];
  o[33] := o[2]*o[21]*o[3]*tau1;
  o[34] := 1/o[33];
  o[35] := o[1]*o[3]*tau1;
  o[36] := 1/o[35];
  o[37] := o[1]*o[3];
  o[38] := 1/o[37];
  o[39] := 1/o[6];
  o[40] := o[1]*o[22]*o[3];
  o[41] := 1/o[40];
  o[42] := 1/o[22];
  o[43] := o[1]*o[2]*o[21]*o[3];
  o[44] := 1/o[43];
  o[45] := 1/o[13];

  g := pi1*(pi1*(pi1*(o[10]*(-0.000031679644845054 + o[2]*(-2.82707979853120e-6
     - 8.5205128120103e-10*o[6])) + pi1*(o[12]*(-2.24252819080000e-6 + (-6.5171222895601e-7
     - 1.43417299379240e-13*o[13])*o[7]) + pi1*(-4.0516996860117e-7*o[14]
     + o[16]*((-1.27343017416410e-9 - 1.74248712306340e-10*o[11])*o[36]
     + o[19]*(-6.8762131295531e-19*o[34] + o[15]*(1.44783078285210e-20*o[
    32] + o[20]*(2.63357816627950e-23*o[30] + pi1*(-1.19476226400710e-23*
    o[28] + pi1*(1.82280945814040e-24*o[26] - 9.3537087292458e-26*o[24]*
    pi1))))))))) + o[8]*(-0.00047184321073267 + o[7]*(-0.000300017807930260
     + (0.000047661393906987 + o[1]*(-4.4141845330846e-6 -
    7.2694996297594e-16*o[9]))*tau1))) + o[5]*(0.000283190801238040 + o[1]
    *(-0.00060706301565874 + o[6]*(-0.0189900682184190 + tau1*(-0.032529748770505
     + (-0.0218417171754140 - 0.000052838357969930*o[1])*tau1))))) + (
    0.146329712131670 + tau1*(-0.84548187169114 + tau1*(-3.7563603672040
     + tau1*(3.3855169168385 + tau1*(-0.95791963387872 + tau1*(
    0.157720385132280 + (-0.0166164171995010 + 0.00081214629983568*tau1)*
    tau1))))))/o[1];
end calc_g;
