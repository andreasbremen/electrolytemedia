within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp.Reduced;
function calc_der_g_tau "Derivative of Gibbs free energy w.r.t. tau"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;

  output Real gtau;
protected
  Real pi1 "Dimensionless pressure";
  Real tau1 "Dimensionless temperature";
  Real[45] o "Vector of auxiliary variables";

  Real pi(unit="1") "Dimensionless pressure";
  Real tau(unit="1") "Dimensionless temperature";
algorithm
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

  gtau := pi1*(o[38]*(-0.00254871721114236 + o[1]*(0.0042494411096112
     + (0.0189900682184190 + (-0.0218417171754140 - 0.000158515073909790*
    o[1])*o[1])*o[6])) + pi1*(o[10]*(0.00141552963219801 + o[2]*(
    0.000047661393906987 + o[1]*(-0.0000132425535992538 -
    1.23581493705910e-14*o[9]))) + pi1*(o[12]*(0.000126718579380216 -
    5.1123076872062e-9*o[37]) + pi1*(o[39]*(0.0000112126409540000 + (
    1.30342445791202e-6 - 1.43417299379240e-12*o[13])*o[7]) + pi1*(
    3.2413597488094e-6*o[5] + o[16]*((1.40077319158051e-8 +
    1.04549227383804e-9*o[11])*o[45] + o[19]*(1.99410180757040e-17*o[44]
     + o[15]*(-4.4882754268415e-19*o[42] + o[20]*(-1.00075970318621e-21*o[
    28] + pi1*(4.6595728296277e-22*o[26] + pi1*(-7.2912378325616e-23*o[24]
     + 3.8350205789908e-24*o[41]*pi1))))))))))) + o[8]*(-0.292659424263340
     + tau1*(0.84548187169114 + o[1]*(3.3855169168385 + tau1*(-1.91583926775744
     + tau1*(0.47316115539684 + (-0.066465668798004 + 0.0040607314991784*
    tau1)*tau1)))));
end calc_der_g_tau;
