within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp.Reduced;
function calc_3der_g_tautaupi
  "Thrid derivative of Gibbs free energy w.r.t. tau,tau,pi"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;

  output Real gtautaupi;
protected
  Real pi1 "Dimensionless pressure";
  Real der_pi1;
  Real tau1 "Dimensionless temperature";
  Real[45] o "Vector of auxiliary variables";
  Real der_o15 "Derivative of auxiliary variable w.r.t. p";
  Real der_o16 "Derivative of auxiliary variable w.r.t. p";
  Real der_o17 "Derivative of auxiliary variable w.r.t. p";
  Real der_o18 "Derivative of auxiliary variable w.r.t. p";
  Real der_o19 "Derivative of auxiliary variable w.r.t. p";
  Real der_o20 "Derivative of auxiliary variable w.r.t. p";

  Real term1;
  Real term2;
  Real term3;
  Real term4;
  Real term5;
  Real term6;
  Real term7;
  Real term8;
  Real term9;
  Real term10;
  Real term11;
  Real term12;
  Real term13;
  Real term14;

  Real der_term2;
  Real der_term4;
  Real der_term6;
  Real der_term8;
  Real der_term10;
  Real der_term11;
  Real der_term12;
  Real der_term13;
  Real der_term14;

  Real pi(unit="1") "Dimensionless pressure";
  Real tau(unit="1") "Dimensionless temperature";
algorithm
  pi := p/IF97.PSTAR1;
  tau := IF97.TSTAR1/T;
  pi1 := 7.1000000000000 - pi;
  der_pi1 :=-1;
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
  der_o15 :=2*pi1*der_pi1;
  o[16] := o[15]*pi1;
  der_o16 :=der_o15*pi1 + o[15]*der_pi1;
  o[17] := o[15]*o[15];
  der_o17 :=2*o[15]*der_o15;
  o[18] := o[17]*o[17];
  der_o18 :=2*o[17]*der_o17;
  o[19] := o[17]*o[18]*pi1;
  der_o19 :=pi1*(o[17]*der_o18 + der_o17*o[18]) + o[17]*o[18]*der_pi1;
  o[20] := o[15]*o[17];
  der_o20 :=der_o15*o[17] + o[15]*der_o17;
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

  term14 :=2.98940751135026e-21*o[41] - (1.61070864317613e-22*pi1)/(o[1]*o[22]*
    o[3]*tau1);
  der_term14 :=-(1.61070864317613e-22*der_pi1)/(o[1]*o[22]*o[3]*tau1);

  term13 :=-1.86382913185108e-20*o[24] + pi1*term14;
  der_term13 :=der_pi1*term14 + pi1*der_term14;

  term12 :=3.9029628424262e-20*o[26] + pi1*term13;
  der_term12 :=der_pi1*term13 + pi1*der_term13;

  term11 :=o[20]*term12 + 1.43624813658928e-17/(o[22]*tau1);
  der_term11 :=der_o20*term12 + o[20]*der_term12;

  term10 :=-5.9823054227112e-16*o[32] + o[15]*term11;
  der_term10 :=der_o15*term11 + o[15]*der_term11;

  term9 :=(-1.68092782989661e-7 - 7.3184459168663e-9*o[11])/(o[2]*o[3]*tau1);

  term8 :=(o[19]*term10 + term9);
  der_term8 :=der_o19*term10 + o[19]*der_term10;

  term7 :=(-0.000067275845724000 + (-3.9102733737361e-6 - 1.29075569441316e-11*
    o[13])*o[7])/(o[1]*o[2]*tau1);

  term6 :=-0.0000291722377392842*o[38] + o[16]*term8;
  der_term6 :=der_o16*term8 + o[16]*der_term8;

  term5 :=(-0.00063359289690108 - 2.55615384360309e-8*o[37])*o[39];

  term4 :=term5 + pi1*(pi1*term6 + term7);
  der_term4 :=der_pi1*(pi1*term6 + term7) + pi1*(der_pi1*term6 + pi1*der_term6);

  term3 :=o[12]*(-0.0056621185287920 + o[6]*(-0.0000264851071985076 - 1.97730389929456e-13
    *o[9]));

  term2 :=term3 + pi1*term4;
  der_term2 :=der_pi1*term4 + pi1*der_term4;

  term1 :=o[36]*(0.0254871721114236 + o[1]*(-0.033995528876889 + (-0.037980136436838
     - 0.00031703014781958*o[2])*o[6]));

  gtautaupi :=der_pi1*(term1 + pi1*term2) + pi1*(der_pi1*term2 + pi1*der_term2);

end calc_3der_g_tautaupi;
