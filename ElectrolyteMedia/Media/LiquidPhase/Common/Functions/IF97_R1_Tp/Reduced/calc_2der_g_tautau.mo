within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp.Reduced;
function calc_2der_g_tautau "Second derivative of Gibbs free energy w.r.t. tau"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;

  output Real gtautau;
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

  gtautau := pi1*(o[36]*(0.0254871721114236 + o[1]*(-0.033995528876889
     + (-0.037980136436838 - 0.00031703014781958*o[2])*o[6])) + pi1*(o[12]
    *(-0.0056621185287920 + o[6]*(-0.0000264851071985076 -
    1.97730389929456e-13*o[9])) + pi1*((-0.00063359289690108 -
    2.55615384360309e-8*o[37])*o[39] + pi1*(pi1*(-0.0000291722377392842*o[
    38] + o[16]*(o[19]*(-5.9823054227112e-16*o[32] + o[15]*(o[20]*(
    3.9029628424262e-20*o[26] + pi1*(-1.86382913185108e-20*o[24] + pi1*(
    2.98940751135026e-21*o[41] - (1.61070864317613e-22*pi1)/(o[1]*o[22]*o[
    3]*tau1)))) + 1.43624813658928e-17/(o[22]*tau1))) + (-1.68092782989661e-7
     - 7.3184459168663e-9*o[11])/(o[2]*o[3]*tau1))) + (-0.000067275845724000
     + (-3.9102733737361e-6 - 1.29075569441316e-11*o[13])*o[7])/(o[1]*o[2]
    *tau1))))) + o[10]*(0.87797827279002 + tau1*(-1.69096374338228 + o[7]
    *(-1.91583926775744 + tau1*(0.94632231079368 + (-0.199397006394012 +
    0.0162429259967136*tau1)*tau1))));
end calc_2der_g_tautau;
