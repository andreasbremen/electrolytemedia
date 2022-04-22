within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp.Reduced;
function calc_2der_g_pipi "Second derivative of Gibbs free energy w.r.t. pi"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;

  output Real gpipi;
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

  gpipi := pi1*(o[10]*(-0.000190077869070324 + o[2]*(-0.0000169624787911872
     - 5.1123076872062e-9*o[6])) + pi1*(o[12]*(-0.0000269103382896000 + (
    -7.8205467474721e-6 - 1.72100759255088e-12*o[13])*o[7]) + pi1*(-8.1033993720234e-6
    *o[14] + o[16]*((-7.1312089753190e-8 - 9.7579278891550e-9*o[11])*o[36]
     + o[19]*(-2.88800951441230e-16*o[34] + o[15]*(7.3260237612316e-18*o[
    32] + o[20]*(2.13846547101895e-20*o[30] + pi1*(-1.03944316968618e-20*
    o[28] + pi1*(1.69521279607057e-21*o[26] - 9.2788790594118e-23*o[24]*
    pi1))))))))) + o[8]*(-0.00094368642146534 + o[7]*(-0.00060003561586052
     + (0.000095322787813974 + o[1]*(-8.8283690661692e-6 -
    1.45389992595188e-15*o[9]))*tau1));
end calc_2der_g_pipi;
