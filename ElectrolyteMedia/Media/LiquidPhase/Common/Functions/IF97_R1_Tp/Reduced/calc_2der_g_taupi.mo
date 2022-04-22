within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp.Reduced;
function calc_2der_g_taupi "Derivative of Gibbs free energy w.r.t. tau and pi"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;

  output Real gtaupi;
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

  gtaupi := o[38]*(0.00254871721114236 + o[1]*(-0.0042494411096112 + (-0.0189900682184190
     + (0.0218417171754140 + 0.000158515073909790*o[1])*o[1])*o[6])) +
    pi1*(o[10]*(-0.00283105926439602 + o[2]*(-0.000095322787813974 + o[1]
    *(0.0000264851071985076 + 2.47162987411820e-14*o[9]))) + pi1*(o[12]*(
    -0.00038015573814065 + 1.53369230616185e-8*o[37]) + pi1*(o[39]*(-0.000044850563816000
     + (-5.2136978316481e-6 + 5.7366919751696e-12*o[13])*o[7]) + pi1*(-0.0000162067987440468
    *o[5] + o[16]*((-1.12061855326441e-7 - 8.3639381907043e-9*o[11])*o[45]
     + o[19]*(-4.1876137958978e-16*o[44] + o[15]*(1.03230334817355e-17*o[
    42] + o[20]*(2.90220313924001e-20*o[28] + pi1*(-1.39787184888831e-20*
    o[26] + pi1*(2.26028372809410e-21*o[24] - 1.22720658527705e-22*o[41]*
    pi1))))))))));
end calc_2der_g_taupi;
