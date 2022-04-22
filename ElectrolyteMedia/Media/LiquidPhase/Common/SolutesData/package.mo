within ElectrolyteMedia.Media.LiquidPhase.Common;
package SolutesData "Properties of fluids from SUPCRTBL by Zimmer(2016) and Johnson(1992)"

  constant DataRecord[:] Acetic = {ACETIC_ACID_aq,ACETATE,OH_n1,H_p1} "Acetic acid, acetate, OH-, H+";
  constant DataRecord[:] NaCl = {Na_p1,Cl_n1} "Na+,Cl-";
  constant DataRecord[:] NaCl_H2O = {Na_p1,Cl_n1,OH_n1,H_p1} "Na+,Cl-,OH-,H+";
  constant DataRecord[:] magnesite_H2O = {Mg_p2,CO3_n2,Na_p1,Cl_n1} "Mg+2,CO3-2,Na+,Cl-";
  constant DataRecord[:] CarbonDioxide = {CO2_aq,HCO3_n1,CO3_n2,OH_n1,H_p1} "CO2,HCO3-,CO3-2,OH-,H+";
  constant DataRecord[:] CarbonDioxide_Acetic = {CO2_aq,HCO3_n1,CO3_n2,ACETIC_ACID_aq,ACETATE,OH_n1,H_p1} "CO2,HCO3-,CO3-2,acetic acid,acetate,OH-,H+";
  constant DataRecord[:] CarbonDioxide_NaCl = {CO2_aq,HCO3_n1,CO3_n2,Na_p1,Cl_n1,OH_n1,H_p1} "CO2,HCO3-,CO3-2,Na+,Cl-,OH-,H+";
  constant DataRecord[:] NaOH = {Na_p1,OH_n1} "Na+,OH-";
  constant DataRecord[:] NaOH_dis = {Na_p1,OH_n1,NaOH_aq} "Na+,OH-,NaOH";
  constant DataRecord[:] NaOH_dis_ = {Na_p1,OH_n1,NaOH_aq,H_p1} "Na+,OH-,NaOH,H+";
  constant DataRecord[:] NaOH_ = {Na_p1,OH_n1,H_p1} "Na+,OH-,H+";
  constant DataRecord[:] Acetic_NaCl = {ACETIC_ACID_aq,ACETATE,Na_p1,Cl_n1,OH_n1,H_p1} "Acetic acid,acetate,Na+,Cl-,OH-,H+";
  constant DataRecord[:] Acetic_NaOH = {ACETIC_ACID_aq,ACETATE,NaOH_aq,Na_p1,Cl_n1,OH_n1,H_p1} "Acetic acid,acetate,NaOH,Na+,Cl-,OH-,H+";
  constant DataRecord[:] CarbonDioxide_Acetic_NaCl = {CO2_aq,HCO3_n1,CO3_n2,ACETIC_ACID_aq,ACETATE,Na_p1,Cl_n1,OH_n1,H_p1} "CO2,HCO3-,CO3-2,acetic acid,acetate,Na+,Cl-,OH-,H+";
  constant DataRecord[:] HydrochloricAcid = {HCl_aq,Na_p1,Cl_n1,OH_n1,H_p1} "HCl,Na+,Cl-,OH-,H+";
  constant DataRecord[:] OxalicAcid = {OXALIC_ACID_aq,H_OXALATE_aq,OXALATE_aq,Na_p1,OH_n1,H_p1} "OXALIC_ACID_aq,H_OXALATE_aq,OXALATE_aq,Na+,OH-,H+";
  constant DataRecord[:] PhosphoricAcid = {H3PO4_aq,H2PO4_n1,HPO4_n2,PO4_n3,Na_p1,OH_n1,H_p1} "H3PO4,aq,H2PO4-,HPO4-2,PO4-3,Na+,OH-,H+";
  constant DataRecord[:] FormicAcid = {FORMIC_ACID_aq,FORMATE_aq,Na_p1,OH_n1,H_p1} "Formic acid,formate,Na+,OH-,H+";
  constant DataRecord[:] KCl_NaCl = {K_p1,Cl_n1,Na_p1,OH_n1,H_p1} "K+,Cl-,Na+,OH-,H+";
  constant DataRecord[:] Seawater = {Na_p1,K_p1,Ca_p2,Mg_p2,Cl_n1,SO4_n2,SiO2_aq,OH_n1,H_p1} "Seawater:Na+,K+,Ca+2,Mg+2,Cl-,SO4-2,SiO2,OH-,H+";
  constant DataRecord[:] Electrolysis = {H_p1, OH_n1,CO_aq,CO2_aq,K_p1,HCO3_n1,CO3_n2,KSO4_n1,SO4_n2,HSO4_n1,KOH_aq,H2_aq,O2_aq}  "Electrolysis:H+,OH-,CO,CO2,K+,HCO3-,CO3-2,KSO4-,SO4-2,HSO4-,KOH,H2,O2";
  constant DataRecord[:] Electrolysis_ = {H_p1, OH_n1,CO_aq,CO2_aq,K_p1,HCO3_n1,CO3_n2,KSO4_n1,SO4_n2,KOH_aq,H2_aq,O2_aq}  "Electrolysis:H+,OH-,CO,CO2,K+,HCO3-,CO3-2,KSO4-,SO4-2,KOH,H2,O2";
  constant DataRecord[:] ReactiveTransport = {CO2_aq,HCO3_n1,CO3_n2,Mg_p2,Ca_p2,Na_p1,Cl_n1,SiO2_aq,OH_n1,H_p1} "ReactiveTransport: CO2,HCO3-,CO3-2,Mg+2,Ca+2,Na+,Cl-,SiO2,OH-,H+";
  constant DataRecord[:] ReactiveTransport_ = {CO2_aq,HCO3_n1,OH_n1,CO3_n2,Mg_p2,Ca_p2,Na_p1,Cl_n1,SiO2_aq,H_p1} "ReactiveTransport: CO2,HCO3-,OH-,CO3-2,Mg+2,Ca+2,Na+,Cl-,SiO2,H+";

  //1
protected 
  constant DataRecord ONE_BUTANAMINE_aq(
    R = Modelica.Constants.R/ONE_BUTANAMINE_aq.MM,
    G_ref = 39371,
    H_ref = -151084,
    S_ref = 197.485,
    a1 = 5.8269e-05,
    a2 = 8251.06,
    a3 = 0.00049989,
    a4 = -150377,
    c1 = 384.5552,
    c2 = 116047,
    w_ref = -156820,
    z = 0,
    a_0 = 3.72,
    MM = 0.073134) "1-BUTANAMINE,aq";

  //2
  constant DataRecord ONE_BUTANOL_aq(
    R = Modelica.Constants.R/ONE_BUTANOL_aq.MM,
    G_ref = -162507,
    H_ref = -336059,
    S_ref = 196.23,
    a1 = 5.6443e-05,
    a2 = 7506.1,
    a3 = 0.00059336,
    a4 = -147298,
    c1 = 397.0022,
    c2 = 125662,
    w_ref = -154890,
    z = 0,
    a_0 = 3.72,
    MM = 0.074119) "1-BUTANOL,aq";

  //3
  constant DataRecord ONE_BUTENE_aq(
    R = Modelica.Constants.R/ONE_BUTENE_aq.MM,
    G_ref = 84977,
    H_ref = -23577,
    S_ref = 181.586,
    a1 = 5.0598e-05,
    a2 = 9099.2,
    a3 = -0.00011711,
    a4 = -153883,
    c1 = 257.69,
    c2 = 767074,
    w_ref = -274970,
    z = 0,
    a_0 = 3.72,
    MM = 0.056103) "1-BUTENE,aq";

  //4
  constant DataRecord ONE_BUTYNE_aq(
    R = Modelica.Constants.R/ONE_BUTYNE_aq.MM,
    G_ref = 209326,
    H_ref = 139746,
    S_ref = 181.586,
    a1 = 5.2774e-05,
    a2 = 9630.94,
    a3 = -0.00013812,
    a4 = -156084,
    c1 = 234.3952,
    c2 = 686109,
    w_ref = -274970,
    z = 0,
    a_0 = 3.72,
    MM = 0.054087) "1-BUTYNE,aq";

  //5
  constant DataRecord ONE_HEPTANAMINE_aq(
    R = Modelica.Constants.R/ONE_HEPTANAMINE_aq.MM,
    G_ref = 70751,
    H_ref = -217526,
    S_ref = 279.91,
    a1 = 8.5214e-05,
    a2 = 13636.53,
    a3 = 0.00054433,
    a4 = -172644,
    c1 = 596.1066,
    c2 = 290842,
    w_ref = -281630,
    z = 0,
    a_0 = 3.72,
    MM = 0.11521) "1-HEPTANAMINE,aq";

  //6
  constant DataRecord ONE_HEPTANOL_aq(
    R = Modelica.Constants.R/ONE_HEPTANOL_aq.MM,
    G_ref = -142653,
    H_ref = -406978,
    S_ref = 301.666,
    a1 = 8.2716e-05,
    a2 = 13136.8,
    a3 = 0.00054028,
    a4 = -170577,
    c1 = 634.217,
    c2 = 323084,
    w_ref = -314590,
    z = 0,
    a_0 = 3.72,
    MM = 0.1162) "1-HEPTANOL,aq";

  //7
  constant DataRecord ONE_HEPTENE_aq(
    R = Modelica.Constants.R/ONE_HEPTENE_aq.MM,
    G_ref = 110667,
    H_ref = -94851,
    S_ref = 265.684,
    a1 = 7.7306e-05,
    a2 = 15620.63,
    a3 = -0.00037346,
    a4 = -180845,
    c1 = 401.9092,
    c2 = 1309128,
    w_ref = -402330,
    z = 0,
    a_0 = 3.72,
    MM = 0.098181) "1-HEPTENE,aq";

  //8
  constant DataRecord ONE_HEPTYNE_aq(
    R = Modelica.Constants.R/ONE_HEPTYNE_aq.MM,
    G_ref = 237358,
    H_ref = 71044,
    S_ref = 266.521,
    a1 = 7.9478e-05,
    a2 = 16150.87,
    a3 = -0.00039432,
    a4 = -183037,
    c1 = 378.498,
    c2 = 1228159,
    w_ref = -403630,
    z = 0,
    a_0 = 3.72,
    MM = 0.096165) "1-HEPTYNE,aq";

  //9
  constant DataRecord ONE_HEXANAMINE_aq(
    R = Modelica.Constants.R/ONE_HEXANAMINE_aq.MM,
    G_ref = 62174,
    H_ref = -193803,
    S_ref = 251.877,
    a1 = 7.6197e-05,
    a2 = 11833.27,
    a3 = 0.00052974,
    a4 = -165189,
    c1 = 532.1642,
    c2 = 237668,
    w_ref = -239200,
    z = 0,
    a_0 = 3.72,
    MM = 0.10119) "1-HEXANAMINE,aq";

  //10
  constant DataRecord ONE_HEXANOL_aq(
    R = Modelica.Constants.R/ONE_HEXANOL_aq.MM,
    G_ref = -155059,
    H_ref = -387815,
    S_ref = 271.123,
    a1 = 7.4409e-05,
    a2 = 11476.25,
    a3 = 0.00052672,
    a4 = -163712,
    c1 = 530.9228,
    c2 = 238798,
    w_ref = -268320,
    z = 0,
    a_0 = 3.72,
    MM = 0.10217) "1-HEXANOL,aq";

  //11
  constant DataRecord ONE_HEXENE_aq(
    R = Modelica.Constants.R/ONE_HEXENE_aq.MM,
    G_ref = 101964,
    H_ref = -71233,
    S_ref = 237.651,
    a1 = 6.8403e-05,
    a2 = 13446.83,
    a3 = -0.00028802,
    a4 = -171858,
    c1 = 353.8363,
    c2 = 1128442,
    w_ref = -359910,
    z = 0,
    a_0 = 3.72,
    MM = 0.084155) "1-HEXENE,aq";

  //12
  constant DataRecord ONE_HEXYNE_aq(
    R = Modelica.Constants.R/ONE_HEXYNE_aq.MM,
    G_ref = 227693,
    H_ref = 93471,
    S_ref = 237.651,
    a1 = 7.0579e-05,
    a2 = 13978.58,
    a3 = -0.00030902,
    a4 = -174059,
    c1 = 330.541,
    c2 = 1047477,
    w_ref = -359910,
    z = 0,
    a_0 = 3.72,
    MM = 0.082139) "1-HEXYNE,aq";

  //13
  constant DataRecord ONE_OCTANAMINE_aq(
    R = Modelica.Constants.R/ONE_OCTANAMINE_aq.MM,
    G_ref = 79329,
    H_ref = -241249,
    S_ref = 307.942,
    a1 = 9.4002e-05,
    a2 = 15392.77,
    a3 = 0.00055887,
    a4 = -179904,
    c1 = 668.7116,
    c2 = 350803,
    w_ref = -324090,
    z = 0,
    a_0 = 3.72,
    MM = 0.12924) "1-OCTANAMINE,aq";

  //14
  constant DataRecord ONE_OCTANOL_aq(
    R = Modelica.Constants.R/ONE_OCTANOL_aq.MM,
    G_ref = -126566,
    H_ref = -431203,
    S_ref = 302.503,
    a1 = 9.174e-05,
    a2 = 14940.77,
    a3 = 0.00055509,
    a4 = -178033,
    c1 = 709.1721,
    c2 = 381916,
    w_ref = -315850,
    z = 0,
    a_0 = 3.72,
    MM = 0.13022) "1-OCTANOL,aq";

  //15
  constant DataRecord ONE_OCTENE_aq(
    R = Modelica.Constants.R/ONE_OCTENE_aq.MM,
    G_ref = 120164,
    H_ref = -117654,
    S_ref = 293.717,
    a1 = 8.6209e-05,
    a2 = 17794.43,
    a3 = -0.00045891,
    a4 = -189832,
    c1 = 449.9821,
    c2 = 1489809,
    w_ref = -444800,
    z = 0,
    a_0 = 3.72,
    MM = 0.11221) "1-OCTENE,aq";

  //16
  constant DataRecord ONE_OCTYNE_aq(
    R = Modelica.Constants.R/ONE_OCTYNE_aq.MM,
    G_ref = 246270,
    H_ref = 47405,
    S_ref = 293.717,
    a1 = 8.8384e-05,
    a2 = 18326.17,
    a3 = -0.00047993,
    a4 = -192029,
    c1 = 426.6877,
    c2 = 1408841,
    w_ref = -444800,
    z = 0,
    a_0 = 3.72,
    MM = 0.11019) "1-OCTYNE,aq";

  //17
  constant DataRecord ONE_PENTANAMINE_aq(
    R = Modelica.Constants.R/ONE_PENTANAMINE_aq.MM,
    G_ref = 53555,
    H_ref = -170080,
    S_ref = 223.844,
    a1 = 6.7237e-05,
    a2 = 10043.61,
    a3 = 0.00051467,
    a4 = -157791,
    c1 = 460.2814,
    c2 = 178272,
    w_ref = -196730,
    z = 0,
    a_0 = 3.72,
    MM = 0.087159) "1-PENTANAMINE,aq";

  //18
  constant DataRecord ONE_PENTANOL_aq(
    R = Modelica.Constants.R/ONE_PENTANOL_aq.MM,
    G_ref = -160958,
    H_ref = -367062,
    S_ref = 223.426,
    a1 = 6.551e-05,
    a2 = 9465.46,
    a3 = 0.00057096,
    a4 = -155398,
    c1 = 468.28,
    c2 = 184493,
    w_ref = -196100,
    z = 0,
    a_0 = 3.72,
    MM = 0.088145) "1-PENTANOL,aq";

  //19
  constant DataRecord ONE_PENTENE_aq(
    R = Modelica.Constants.R/ONE_PENTENE_aq.MM,
    G_ref = 94014,
    H_ref = -46861,
    S_ref = 209.618,
    a1 = 5.9501e-05,
    a2 = 11272.99,
    a3 = -0.00020256,
    a4 = -162871,
    c1 = 305.763,
    c2 = 947760,
    w_ref = -317440,
    z = 0,
    a_0 = 3.72,
    MM = 0.070129) "1-PENTENE,aq";

  //20
  constant DataRecord ONE_PENTYNE_aq(
    R = Modelica.Constants.R/ONE_PENTYNE_aq.MM,
    G_ref = 218237,
    H_ref = 116315,
    S_ref = 209.618,
    a1 = 6.1676e-05,
    a2 = 11804.74,
    a3 = -0.00022356,
    a4 = -165071,
    c1 = 282.469,
    c2 = 866791,
    w_ref = -317440,
    z = 0,
    a_0 = 3.72,
    MM = 0.068113) "1-PENTYNE,aq";

  //21
  constant DataRecord ONE_PROPANAMINE_aq(
    R = Modelica.Constants.R/ONE_PROPANAMINE_aq.MM,
    G_ref = 29330,
    H_ref = -128365,
    S_ref = 170.707,
    a1 = 4.9446e-05,
    a2 = 6485.74,
    a3 = 0.00048582,
    a4 = -143080,
    c1 = 305.2784,
    c2 = 50995,
    w_ref = -116270,
    z = 0,
    a_0 = 3.72,
    MM = 0.059108) "1-PROPANAMINE,aq";

  //22
  constant DataRecord ONE_PROPANOL_aq(
    R = Modelica.Constants.R/ONE_PROPANOL_aq.MM,
    G_ref = -175351,
    H_ref = -315139,
    S_ref = 173.218,
    a1 = 4.7457e-05,
    a2 = 5993.66,
    a3 = 0.0005066,
    a4 = -141047,
    c1 = 327.6666,
    c2 = 68814,
    w_ref = -120040,
    z = 0,
    a_0 = 3.72,
    MM = 0.060093) "1-PROPANOL,aq";

  //23
  constant DataRecord ONE_PROPENE_aq(
    R = Modelica.Constants.R/ONE_PROPENE_aq.MM,
    G_ref = 74935,
    H_ref = -1213,
    S_ref = 153.553,
    a1 = 4.1696e-05,
    a2 = 6925.27,
    a3 = -3.1648e-05,
    a4 = -144900,
    c1 = 209.6171,
    c2 = 586392,
    w_ref = -232550,
    z = 0,
    a_0 = 3.72,
    MM = 0.042077) "1-PROPENE,aq";

  //24
  constant DataRecord ONE_PROPYNE_aq(
    R = Modelica.Constants.R/ONE_PROPYNE_aq.MM,
    G_ref = 200330,
    H_ref = 163050,
    S_ref = 153.553,
    a1 = 4.3871e-05,
    a2 = 7457.02,
    a3 = -5.2668e-05,
    a4 = -147097,
    c1 = 186.3223,
    c2 = 505423,
    w_ref = -232550,
    z = 0,
    a_0 = 3.72,
    MM = 0.040062) "1-PROPYNE,aq";

  //25
  constant DataRecord TWO_BUTANONE_aq(
    R = Modelica.Constants.R/TWO_BUTANONE_aq.MM,
    G_ref = -153678,
    H_ref = -284010,
    S_ref = 210.455,
    a1 = 5.4034e-05,
    a2 = 7404.3,
    a3 = 0.00049304,
    a4 = -146879,
    c1 = 308.3955,
    c2 = 57781,
    w_ref = -176440,
    z = 0,
    a_0 = 3.72,
    MM = 0.072103) "2-BUTANONE,aq";

  //26
  constant DataRecord TWO_HEPTANONE_aq(
    R = Modelica.Constants.R/TWO_HEPTANONE_aq.MM,
    G_ref = -127319,
    H_ref = -355180,
    S_ref = 293.298,
    a1 = 8.0508e-05,
    a2 = 12693.92,
    a3 = 0.00053708,
    a4 = -168745,
    c1 = 526.3853,
    c2 = 237668,
    w_ref = -301920,
    z = 0,
    a_0 = 3.72,
    MM = 0.11418) "2-HEPTANONE,aq";

  //27
  constant DataRecord TWO_HEXANONE_aq(
    R = Modelica.Constants.R/TWO_HEXANONE_aq.MM,
    G_ref = -135896,
    H_ref = -331456,
    S_ref = 265.266,
    a1 = 7.1605e-05,
    a2 = 10915.47,
    a3 = 0.00052219,
    a4 = -161394,
    c1 = 453.7799,
    c2 = 177707,
    w_ref = -259450,
    z = 0,
    a_0 = 3.72,
    MM = 0.10015) "2-HEXANONE,aq";

  //28
  constant DataRecord TWO_HYDROXYBUTANOATE(
    R = Modelica.Constants.R/TWO_HYDROXYBUTANOATE.MM,
    G_ref = -504381,
    H_ref = -710485,
    S_ref = 161.502,
    a1 = 4.9862e-05,
    a2 = 8481.8,
    a3 = 1.013e-06,
    a4 = -151331,
    c1 = 270.6642,
    c2 = -147503,
    w_ref = 437190,
    z = -1,
    a_0 = 3.72,
    MM = 0.1031) "2-HYDROXYBUTANOATE";

  //29
  constant DataRecord TWO_HYDROXYBUTANOIC(
    R = Modelica.Constants.R/TWO_HYDROXYBUTANOIC.MM,
    G_ref = -526138,
    H_ref = -709899,
    S_ref = 236.396,
    a1 = 5.6007e-05,
    a2 = 9328.27,
    a3 = 0.00010804,
    a4 = -154833,
    c1 = 327.9122,
    c2 = 155791,
    w_ref = -89040,
    z = 0,
    a_0 = 3.72,
    MM = 0.1041) "2-HYDROXYBUTANOIC";

  //30
  constant DataRecord TWO_HYDROXYDECANOATE(
    R = Modelica.Constants.R/TWO_HYDROXYDECANOATE.MM,
    G_ref = -453462,
    H_ref = -853243,
    S_ref = 329.699,
    a1 = 0.00010304,
    a2 = 20540.3,
    a3 = -0.00027416,
    a4 = -201183,
    c1 = 766.277,
    c2 = -82998,
    w_ref = 182260,
    z = -1,
    a_0 = 3.72,
    MM = 0.18725) "2-HYDROXYDECANOATE";

  //31
  constant DataRecord TWO_HYDROXYDECANOIC(
    R = Modelica.Constants.R/TWO_HYDROXYDECANOIC.MM,
    G_ref = -474800,
    H_ref = -852239,
    S_ref = 404.593,
    a1 = 0.00011023,
    a2 = 22168.84,
    a3 = -0.00031093,
    a4 = -207916,
    c1 = 763.5382,
    c2 = 681611,
    w_ref = 22470,
    z = 0,
    a_0 = 3.72,
    MM = 0.18826) "2-HYDROXYDECANOIC";

  //32
  constant DataRecord TWO_HYDROXYHEPTANOIC(
    R = Modelica.Constants.R/TWO_HYDROXYHEPTANOIC.MM,
    G_ref = -500406,
    H_ref = -781069,
    S_ref = 320.494,
    a1 = 8.2905e-05,
    a2 = 15977.36,
    a3 = -0.00017078,
    a4 = -181481,
    c1 = 545.7254,
    c2 = 418701,
    w_ref = -33260,
    z = 0,
    a_0 = 3.72,
    MM = 0.14618) "2-HYDROXYHEPTANOIC";

  //33
  constant DataRecord TWO_HYDROXYHEXANOATE(
    R = Modelica.Constants.R/TWO_HYDROXYHEXANOATE.MM,
    G_ref = -487645,
    H_ref = -758308,
    S_ref = 217.568,
    a1 = 6.7735e-05,
    a2 = 12731.03,
    a3 = -0.00014132,
    a4 = -168900,
    c1 = 435.8632,
    c2 = -126001,
    w_ref = 352170,
    z = -1,
    a_0 = 3.72,
    MM = 0.13115) "2-HYDROXYHEXANOATE";

  //34
  constant DataRecord TWO_HYDROXYHEXANOIC(
    R = Modelica.Constants.R/TWO_HYDROXYHEXANOIC.MM,
    G_ref = -508984,
    H_ref = -757346,
    S_ref = 292.462,
    a1 = 7.4346e-05,
    a2 = 13974.31,
    a3 = -0.00011069,
    a4 = -174038,
    c1 = 473.1213,
    c2 = 331063,
    w_ref = -51880,
    z = 0,
    a_0 = 3.72,
    MM = 0.13215) "2-HYDROXYHEXANOIC";

  //35
  constant DataRecord TWO_HYDROXYNONANOATE(
    R = Modelica.Constants.R/TWO_HYDROXYNONANOATE.MM,
    G_ref = -461997,
    H_ref = -829478,
    S_ref = 301.666,
    a1 = 9.414e-05,
    a2 = 18522.57,
    a3 = -0.00022831,
    a4 = -192841,
    c1 = 683.6744,
    c2 = -93747,
    w_ref = 224760,
    z = -1,
    a_0 = 3.72,
    MM = 0.17322) "2-HYDROXYNONANOATE";

  //36
  constant DataRecord TWO_HYDROXYNONANOIC(
    R = Modelica.Constants.R/TWO_HYDROXYNONANOIC.MM,
    G_ref = -483336,
    H_ref = -828516,
    S_ref = 376.56,
    a1 = 0.00010112,
    a2 = 20105.92,
    a3 = -0.00026443,
    a4 = -199389,
    c1 = 690.934,
    c2 = 593973,
    w_ref = 3890,
    z = 0,
    a_0 = 3.72,
    MM = 0.17423) "2-HYDROXYNONANOIC";

  //37
  constant DataRecord TWO_HYDROXYOCTANOATE(
    R = Modelica.Constants.R/TWO_HYDROXYOCTANOATE.MM,
    G_ref = -470491,
    H_ref = -805713,
    S_ref = 273.634,
    a1 = 8.5237e-05,
    a2 = 16504.79,
    a3 = -0.00018245,
    a4 = -184502,
    c1 = 601.0701,
    c2 = -104500,
    w_ref = 267230,
    z = -1,
    a_0 = 3.72,
    MM = 0.1592) "2-HYDROXYOCTANOATE";

  //38
  constant DataRecord TWO_HYDROXYOCTANOIC(
    R = Modelica.Constants.R/TWO_HYDROXYOCTANOIC.MM,
    G_ref = -491829,
    H_ref = -804792,
    S_ref = 348.527,
    a1 = 9.2014e-05,
    a2 = 18040.28,
    a3 = -0.00021727,
    a4 = -190849,
    c1 = 618.3295,
    c2 = 506339,
    w_ref = -14690,
    z = 0,
    a_0 = 3.72,
    MM = 0.16021) "2-HYDROXYOCTANOIC";

  //39
  constant DataRecord TWO_HYDROXYPENTANOIC(
    R = Modelica.Constants.R/TWO_HYDROXYPENTANOIC.MM,
    G_ref = -517561,
    H_ref = -733623,
    S_ref = 264.429,
    a1 = 6.467e-05,
    a2 = 11312.32,
    a3 = 5.8174e-05,
    a4 = -163034,
    c1 = 400.5163,
    c2 = 243429,
    w_ref = -70460,
    z = 0,
    a_0 = 3.72,
    MM = 0.11813) "2-HYDROXYPENTANOIC";

  //40
  constant DataRecord TWO_OCTANONE_aq(
    R = Modelica.Constants.R/TWO_OCTANONE_aq.MM,
    G_ref = -118742,
    H_ref = -378903,
    S_ref = 321.331,
    a1 = 8.9411e-05,
    a2 = 14474.97,
    a3 = 0.00055132,
    a4 = -176109,
    c1 = 598.9898,
    c2 = 297629,
    w_ref = -344390,
    z = 0,
    a_0 = 3.72,
    MM = 0.12821) "2-OCTANONE,aq";

  //41
  constant DataRecord TWO_PENTANONE_aq(
    R = Modelica.Constants.R/TWO_PENTANONE_aq.MM,
    G_ref = -143888,
    H_ref = -307357,
    S_ref = 235.559,
    a1 = 6.2769e-05,
    a2 = 9148.44,
    a3 = 0.00050782,
    a4 = -154088,
    c1 = 381.4084,
    c2 = 117746,
    w_ref = -214470,
    z = 0,
    a_0 = 3.72,
    MM = 0.086129) "2-PENTANONE,aq";

  //42
  constant DataRecord A_AMINOBUTYRIC_aq(
    R = Modelica.Constants.R/A_AMINOBUTYRIC_aq.MM,
    G_ref = -364510,
    H_ref = -578145,
    S_ref = 195.393,
    a1 = 5.0665e-05,
    a2 = 6730.47,
    a3 = 0.00048756,
    a4 = -144093,
    c1 = 211.9635,
    c2 = -19426,
    w_ref = -153640,
    z = 0,
    a_0 = 3.72,
    MM = 0.10312) "A-AMINOBUTYRIC,aq";

  //43
  constant DataRecord ACETAMIDE_aq(
    R = Modelica.Constants.R/ACETAMIDE_aq.MM,
    G_ref = -212798,
    H_ref = -323381,
    S_ref = 165.268,
    a1 = 3.8882e-06,
    a2 = 5999.31,
    a3 = 5.6229e-05,
    a4 = -141072,
    c1 = 157.3711,
    c2 = -49652,
    w_ref = -136190,
    z = 0,
    a_0 = 3.72,
    MM = 0.059066) "ACETAMIDE,aq";

  //44
  constant DataRecord ACETATE(
    R = Modelica.Constants.R/ACETATE.MM,
    G_ref = -369322,
    H_ref = -486013,
    S_ref = 86.19,
    a1 = 3.2437e-05,
    a2 = 3639.91,
    a3 = 0.00031725,
    a4 = -131315,
    c1 = 110.0392,
    c2 = -161502,
    w_ref = 551530,
    z = -1,
    a_0 = 3.72,
    MM = 0.059044) "ACETATE";

  //45
  constant DataRecord ACETIC_ACID_aq(
    R = Modelica.Constants.R/ACETIC_ACID_aq.MM,
    G_ref = -396476,
    H_ref = -485762,
    S_ref = 178.657,
    a1 = 4.8617e-05,
    a2 = 2183.21,
    a3 = 0.00010497,
    a4 = -125294,
    c1 = 176.046,
    c2 = -64505,
    w_ref = -62760,
    z = 0,
    a_0 = 3.72,
    MM = 0.060052) "ACETIC-ACID,aq";

  //46
  constant DataRecord ACETONE_aq(
    R = Modelica.Constants.R/ACETONE_aq.MM,
    G_ref = -161084,
    H_ref = -258236,
    S_ref = 185.77,
    a1 = 4.5229e-05,
    a2 = 5644.72,
    a3 = 0.00047844,
    a4 = -139603,
    c1 = 229.9091,
    c2 = -6418,
    w_ref = -139080,
    z = 0,
    a_0 = 3.72,
    MM = 0.058077) "ACETONE,aq";

  //47
  constant DataRecord ADIPATE_aq(
    R = Modelica.Constants.R/ADIPATE_aq.MM,
    G_ref = -666888,
    H_ref = -953032,
    S_ref = 137.654,
    a1 = 6.1868e-05,
    a2 = 11208.68,
    a3 = -6.2245e-05,
    a4 = -162607,
    c1 = 42.7726,
    c2 = -183807,
    w_ref = 1135160,
    z = -2,
    a_0 = 3.72,
    MM = 0.14412) "ADIPATE,aq";

  //48
  constant DataRecord ADIPIC_ACID_aq(
    R = Modelica.Constants.R/ADIPIC_ACID_aq.MM,
    G_ref = -722953,
    H_ref = -961274,
    S_ref = 340.159,
    a1 = 7.3244e-05,
    a2 = 13787.33,
    a3 = -0.0001209,
    a4 = -173268,
    c1 = 308.1667,
    c2 = 123549,
    w_ref = -20250,
    z = 0,
    a_0 = 3.72,
    MM = 0.14614) "ADIPIC-ACID,aq";

  //49
  constant DataRecord ALANINE_aq(
    R = Modelica.Constants.R/ALANINE_aq.MM,
    G_ref = -371539,
    H_ref = -552832,
    S_ref = 167.36,
    a1 = 4.1619e-05,
    a2 = 4921.6,
    a3 = 0.00047289,
    a4 = -136616,
    c1 = 146.268,
    c2 = -74270,
    w_ref = -111210,
    z = 0,
    a_0 = 3.72,
    MM = 0.089092) "ALANINE,aq";

  //50
  constant DataRecord ALANATE_aq(
    R = Modelica.Constants.R/ALANATE_aq.MM,
    G_ref = -312377,
    H_ref = -505887,
    S_ref = 126.566,
    a1 = 4.4468e-05,
    a2 = 7265.43,
    a3 = 2.7217e-05,
    a4 = -146306,
    c1 = 53.074,
    c2 = -175138,
    w_ref = 489530,
    z = -1,
    a_0 = 3.72,
    MM = 0.088084) "ALANATE,aq";

  //51
  constant DataRecord ALANYLGLYCINE_aq(
    R = Modelica.Constants.R/ALANYLGLYCINE_aq.MM,
    G_ref = -488398,
    H_ref = -778684,
    S_ref = 212.129,
    a1 = 6.1297e-05,
    a2 = 8855.02,
    a3 = 0.00050522,
    a4 = -152875,
    c1 = 229.0213,
    c2 = 33384,
    w_ref = -178990,
    z = 0,
    a_0 = 3.72,
    MM = 0.14614) "ALANYLGLYCINE,aq";

  //52
  constant DataRecord ASPARAGINE_aq(
    R = Modelica.Constants.R/ASPARAGINE_aq.MM,
    G_ref = -538272,
    H_ref = -780985,
    S_ref = 230.957,
    a1 = 5.0872e-05,
    a2 = 6770.21,
    a3 = 0.00048834,
    a4 = -144256,
    c1 = 123.2673,
    c2 = -85048,
    w_ref = -207480,
    z = 0,
    a_0 = 3.72,
    MM = 0.13212) "ASPARAGINE,aq";

  //53
  constant DataRecord ASPARTIC_ACID_aq(
    R = Modelica.Constants.R/ASPARTIC_ACID_aq.MM,
    G_ref = -721322,
    H_ref = -947132,
    S_ref = 229.283,
    a1 = 4.7795e-05,
    a2 = 6155.67,
    a3 = 0.00048313,
    a4 = -141716,
    c1 = 125.3058,
    c2 = -83634,
    w_ref = -204970,
    z = 0,
    a_0 = 3.72,
    MM = 0.1331) "ASPARTIC-ACID,aq";

  //54
  constant DataRecord AZELAIC_ACID_aq(
    R = Modelica.Constants.R/AZELAIC_ACID_aq.MM,
    G_ref = -684753,
    H_ref = -1007089,
    S_ref = 425.094,
    a1 = 0.00010014,
    a2 = 19884.59,
    a3 = -0.00025954,
    a4 = -198472,
    c1 = 526.366,
    c2 = 386869,
    w_ref = 36070,
    z = 0,
    a_0 = 3.72,
    MM = 0.18822) "AZELAIC-ACID,aq";

  //55
  constant DataRecord AZELATE_aq(
    R = Modelica.Constants.R/AZELATE_aq.MM,
    G_ref = -628144,
    H_ref = -1011105,
    S_ref = 221.752,
    a1 = 8.8545e-05,
    a2 = 17254.48,
    a3 = -0.00019951,
    a4 = -187598,
    c1 = 291.0813,
    c2 = -151507,
    w_ref = 1008760,
    z = -2,
    a_0 = 3.72,
    MM = 0.1862) "AZELATE,aq";

  //56
  constant DataRecord Ag_CO3_n1(
    R = Modelica.Constants.R/Ag_CO3_n1.MM,
    G_ref = -466223,
    H_ref = -596178,
    S_ref = 0,
    a1 = 9.4508e-06,
    a2 = -946.92,
    a3 = 0.00027751,
    a4 = -112357,
    c1 = 63.6583,
    c2 = -213899,
    w_ref = 682410,
    z = -1,
    a_0 = 3.72,
    MM = 0.16788) "Ag(CO3)-";

  //57
  constant DataRecord Ag_CO32_n3(
    R = Modelica.Constants.R/Ag_CO32_n3.MM,
    G_ref = -991148,
    H_ref = -1272773,
    S_ref = -79.496,
    a1 = 1.5534e-05,
    a2 = 538.56,
    a3 = 0.00021912,
    a4 = -118499,
    c1 = 152.9624,
    c2 = -368163,
    w_ref = 2133510,
    z = -3,
    a_0 = 3.72,
    MM = 0.22789) "Ag(CO3)2-3";

  //58
  constant DataRecord Ag_p1(
    R = Modelica.Constants.R/Ag_p1.MM,
    G_ref = 77099,
    H_ref = 105751,
    S_ref = 73.387,
    a1 = 7.232e-06,
    a2 = -1489.84,
    a3 = 0.00029914,
    a4 = -110115,
    c1 = 53.4975,
    c2 = -59639,
    w_ref = 90370,
    z = 1,
    a_0 = 2.5,
    MM = 0.10787) "Ag+";

  //59
  constant DataRecord Ag_p2(
    R = Modelica.Constants.R/Ag_p2.MM,
    G_ref = 269031,
    H_ref = 268613,
    S_ref = -87.864,
    a1 = -1.7987e-06,
    a2 = -3694.89,
    a3 = 0.00038581,
    a4 = -100998,
    c1 = 63.6905,
    c2 = -172138,
    w_ref = 552330,
    z = 2,
    a_0 = 3.72,
    MM = 0.10787) "Ag+2";

  //60
  constant DataRecord AgCl_aq(
    R = Modelica.Constants.R/AgCl_aq.MM,
    G_ref = -73011,
    H_ref = -76442,
    S_ref = 142.674,
    a1 = 2.1794e-05,
    a2 = 2066.85,
    a3 = 0.00015906,
    a4 = -124817,
    c1 = 41.0735,
    c2 = -69864,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14332) "AgCl,aq";

  //61
  constant DataRecord AgCl2_n1(
    R = Modelica.Constants.R/AgCl2_n1.MM,
    G_ref = -215727,
    H_ref = -255768,
    S_ref = 196.648,
    a1 = 3.981e-05,
    a2 = 6466.12,
    a3 = -1.3857e-05,
    a4 = -143001,
    c1 = 80.27,
    c2 = -60488,
    w_ref = 383630,
    z = -1,
    a_0 = 3.72,
    MM = 0.17877) "AgCl2-";

  //62
  constant DataRecord AgCl3_n2(
    R = Modelica.Constants.R/AgCl3_n2.MM,
    G_ref = -346059,
    H_ref = -441885,
    S_ref = 186.188,
    a1 = 6.0685e-05,
    a2 = 11563.03,
    a3 = -0.00021419,
    a4 = -164071,
    c1 = 149.929,
    c2 = -35857,
    w_ref = 1062820,
    z = -2,
    a_0 = 3.72,
    MM = 0.21422) "AgCl3-2";

  //63
  constant DataRecord AgCl4_n3(
    R = Modelica.Constants.R/AgCl4_n3.MM,
    G_ref = -469780,
    H_ref = -630177,
    S_ref = 146.44,
    a1 = 8.3841e-05,
    a2 = 17217.33,
    a3 = -0.00043643,
    a4 = -187447,
    c1 = 230.2196,
    c2 = 10167,
    w_ref = 1790580,
    z = -3,
    a_0 = 3.72,
    MM = 0.24967) "AgCl4-3";

  //64
  constant DataRecord AgOH_aq(
    R = Modelica.Constants.R/AgOH_aq.MM,
    G_ref = -91630,
    H_ref = -133051,
    S_ref = 71.965,
    a1 = 9.5751e-06,
    a2 = -916.8,
    a3 = 0.00027641,
    a4 = -112478,
    c1 = -1.3477,
    c2 = -217304,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.12488) "AgOH,aq";

  //65
  constant DataRecord AgO_n1(
    R = Modelica.Constants.R/AgO_n1.MM,
    G_ref = -23012,
    H_ref = -68199,
    S_ref = 58.158,
    a1 = 6.1254e-06,
    a2 = -1757.03,
    a3 = 0.00030894,
    a4 = -109006,
    c1 = 3.6761,
    c2 = -393727,
    w_ref = 592910,
    z = -1,
    a_0 = 3.72,
    MM = 0.12387) "AgO-";

  //66
  constant DataRecord AgNO3_aq(
    R = Modelica.Constants.R/AgNO3_aq.MM,
    G_ref = -32376,
    H_ref = -97571,
    S_ref = 205.016,
    a1 = 2.8309e-05,
    a2 = 3657.78,
    a3 = 9.6521e-05,
    a4 = -131394,
    c1 = 101.0373,
    c2 = 139628,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.16988) "AgNO3,aq";

  //67
  constant DataRecord AlOH_p2(
    R = Modelica.Constants.R/AlOH_p2.MM,
    G_ref = -696322,
    H_ref = -808822,
    S_ref = -181.125,
    a1 = -1.8962e-06,
    a2 = -3718.66,
    a3 = 0.00038674,
    a4 = -100901,
    c1 = 64.4884,
    c2 = -203418,
    w_ref = 665130,
    z = 2,
    a_0 = 3.72,
    MM = 0.043988) "AlOH+2";

  //68
  constant DataRecord Al_p3(
    R = Modelica.Constants.R/Al_p3.MM,
    G_ref = -487478,
    H_ref = -538761,
    S_ref = -339.753,
    a1 = -1.4219e-05,
    a2 = -6727.41,
    a3 = 0.00050501,
    a4 = -88462,
    c1 = 60.373,
    c2 = -370380,
    w_ref = 1146540,
    z = 3,
    a_0 = 9,
    MM = 0.02698) "Al+3";

  //69
  constant DataRecord AlH3SiO4_p2(
    R = Modelica.Constants.R/AlH3SiO4_p2.MM,
    G_ref = -1781798,
    H_ref = -1959668,
    S_ref = 56.777,
    a1 = 6.694e-07,
    a2 = -3025.03,
    a3 = 0.00036024,
    a4 = -103763,
    c1 = 155.1009,
    c2 = -2077774,
    w_ref = 368190,
    z = 2,
    a_0 = 3.72,
    MM = 0.12208) "AlH3SiO4+2";

  //70
  constant DataRecord Al_OH3_aq(
    R = Modelica.Constants.R/Al_OH3_aq.MM,
    G_ref = -1101735,
    H_ref = -1242656,
    S_ref = 59.35,
    a1 = 2.2855e-05,
    a2 = 2324.63,
    a3 = 0.00014921,
    a4 = -125884,
    c1 = 83.793,
    c2 = 74597,
    w_ref = 0,
    z = 0,
    a_0 = 3.72,
    MM = 0.078004) "Al(OH)3,aq";

  //71
  constant DataRecord Al_OH4_n1(
    R = Modelica.Constants.R/Al_OH4_n1.MM,
    G_ref = -1305772,
    H_ref = -1483575,
    S_ref = 103.546,
    a1 = 3.5538e-05,
    a2 = 5421.46,
    a3 = 2.7489e-05,
    a4 = -138687,
    c1 = 233.1597,
    c2 = -477173,
    w_ref = 435260,
    z = -1,
    a_0 = 3.72,
    MM = 0.095012) "Al(OH)4-";

  //72
  constant DataRecord AlO_p1(
    R = Modelica.Constants.R/AlO_p1.MM,
    G_ref = -661859,
    H_ref = -715046,
    S_ref = -112.968,
    a1 = 9.0814e-06,
    a2 = -1038.09,
    a3 = 0.00028134,
    a4 = -111976,
    c1 = -10.8713,
    c2 = -382648,
    w_ref = 400410,
    z = 1,
    a_0 = 3.72,
    MM = 0.04298) "AlO+";

  //73
  constant DataRecord Al_OH2_p1(
    R = Modelica.Constants.R/Al_OH2_p1.MM,
    G_ref = -899506,
    H_ref = -1016277,
    S_ref = -27.547,
    a1 = 1.0437e-05,
    a2 = -707.47,
    a3 = 0.00026839,
    a4 = -113349,
    c1 = 70.0565,
    c2 = -43786,
    w_ref = 222760,
    z = 1,
    a_0 = 3.72,
    MM = 0.060996) "Al(OH)2+";

  //74
  constant DataRecord HAlO2_aq(
    R = Modelica.Constants.R/HAlO2_aq.MM,
    G_ref = -869017,
    H_ref = -951860,
    S_ref = 20.92,
    a1 = 1.4785e-05,
    a2 = 355.01,
    a3 = 0.00022649,
    a4 = -117738,
    c1 = -97.9596,
    c2 = -553104,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.059988) "HAlO2,aq";

  //75
  constant DataRecord AlO2_n1(
    R = Modelica.Constants.R/AlO2_n1.MM,
    G_ref = -831332,
    H_ref = -929371,
    S_ref = -30.208,
    a1 = 1.5573e-05,
    a2 = 1671.68,
    a3 = -6.6438e-05,
    a4 = -123181,
    c1 = 63.7604,
    c2 = -228384,
    w_ref = 728770,
    z = -1,
    a_0 = 3.72,
    MM = 0.05898) "AlO2-";

  //76
  constant DataRecord Au_p1(
    R = Modelica.Constants.R/Au_p1.MM,
    G_ref = 163176,
    H_ref = 200414,
    S_ref = 107.11,
    a1 = 1.4774e-05,
    a2 = 352.63,
    a3 = 0.00022652,
    a4 = -117725,
    c1 = 31.4172,
    c2 = -129520,
    w_ref = 68950,
    z = 1,
    a_0 = 3.72,
    MM = 0.19697) "Au+";

  //77
  constant DataRecord Au_p3(
    R = Modelica.Constants.R/Au_p3.MM,
    G_ref = 433462,
    H_ref = 409195,
    S_ref = -229.283,
    a1 = -7.1827e-06,
    a2 = -5006.32,
    a3 = 0.00043661,
    a4 = -95575,
    c1 = 98.6483,
    c2 = -196849,
    w_ref = 1008970,
    z = 3,
    a_0 = 3.72,
    MM = 0.19697) "Au+3";

  //78
  constant DataRecord B_OH3_aq(
    R = Modelica.Constants.R/B_OH3_aq.MM,
    G_ref = -968763,
    H_ref = -1074535,
    S_ref = 154.808,
    a1 = 2.9557e-05,
    a2 = 3704.81,
    a3 = 0.00014997,
    a4 = -131591,
    c1 = 47.5169,
    c2 = -24694,
    w_ref = -83680,
    z = 0,
    a_0 = 3.72,
    MM = 0.061834) "B(OH)3,aq";

  //79
  constant DataRecord BENZENE_aq(
    R = Modelica.Constants.R/BENZENE_aq.MM,
    G_ref = 133888,
    H_ref = 51170,
    S_ref = 148.532,
    a1 = 5.4912e-05,
    a2 = 7579.02,
    a3 = 0.00049465,
    a4 = -147599,
    c1 = 338.3291,
    c2 = 74471,
    w_ref = -82680,
    z = 0,
    a_0 = 3.72,
    MM = 0.078107) "BENZENE,aq";

  //80
  constant DataRecord BENZOATE_aq(
    R = Modelica.Constants.R/BENZOATE_aq.MM,
    G_ref = -210874,
    H_ref = -355598,
    S_ref = 151.461,
    a1 = 5.8742e-05,
    a2 = 10466.4,
    a3 = -3.7597e-05,
    a4 = -159536,
    c1 = 286.3906,
    c2 = -145720,
    w_ref = 452420,
    z = -1,
    a_0 = 3.72,
    MM = 0.12111) "BENZOATE,aq";

  //81
  constant DataRecord BENZOIC_ACID_aq(
    R = Modelica.Constants.R/BENZOIC_ACID_aq.MM,
    G_ref = -234848,
    H_ref = -355933,
    S_ref = 230.538,
    a1 = 6.362e-05,
    a2 = 11606.08,
    a3 = -7.1333e-05,
    a4 = -164247,
    c1 = 331.2326,
    c2 = 160339,
    w_ref = -92930,
    z = 0,
    a_0 = 3.72,
    MM = 0.12212) "BENZOIC-ACID,aq";

  //82
  constant DataRecord BF4_n1(
    R = Modelica.Constants.R/BF4_n1.MM,
    G_ref = -1486994,
    H_ref = -1574858,
    S_ref = 179.912,
    a1 = 3.3387e-05,
    a2 = 4896.41,
    a3 = 4.8133e-05,
    a4 = -136516,
    c1 = 49.7649,
    c2 = -174695,
    w_ref = 409150,
    z = -1,
    a_0 = 3.72,
    MM = 0.08681) "BF4-";

  //83
  constant DataRecord BO2_n1(
    R = Modelica.Constants.R/BO2_n1.MM,
    G_ref = -678812,
    H_ref = -772366,
    S_ref = -37.238,
    a1 = -9.3839e-06,
    a2 = -2596.8,
    a3 = -0.0002645,
    a4 = -105537,
    c1 = -6.9124,
    c2 = -476403,
    w_ref = 736170,
    z = -1,
    a_0 = 3.72,
    MM = 0.04281) "BO2-";

  //84
  constant DataRecord BUTANE_aq(
    R = Modelica.Constants.R/BUTANE_aq.MM,
    G_ref = 151,
    H_ref = -151586,
    S_ref = 167.444,
    a1 = 5.3934e-05,
    a2 = 9914.41,
    a3 = -0.0001493,
    a4 = -157256,
    c1 = 330.7741,
    c2 = 1014235,
    w_ref = -253590,
    z = 0,
    a_0 = 3.72,
    MM = 0.058119) "BUTANE,aq";

  //85
  constant DataRecord BUTANOIC_ACID_aq(
    R = Modelica.Constants.R/BUTANOIC_ACID_aq.MM,
    G_ref = -381623,
    H_ref = -535343,
    S_ref = 234.722,
    a1 = 5.5522e-05,
    a2 = 9769.1,
    a3 = -2.9225e-05,
    a4 = -156653,
    c1 = 303.0647,
    c2 = 125202,
    w_ref = -90170,
    z = 0,
    a_0 = 3.72,
    MM = 0.088103) "BUTANOIC-ACID,aq";

  //86
  constant DataRecord Ba_CO3_aq(
    R = Modelica.Constants.R/Ba_CO3_aq.MM,
    G_ref = -1103865,
    H_ref = -1195996,
    S_ref = 66.944,
    a1 = 1.3096e-06,
    a2 = -2936.33,
    a3 = 0.00035602,
    a4 = -104140,
    c1 = -60.9944,
    c2 = -423546,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.19735) "Ba(CO3),aq";

  //87
  constant DataRecord Ba_HCO3_p1(
    R = Modelica.Constants.R/Ba_HCO3_p1.MM,
    G_ref = -1153529,
    H_ref = -1231560,
    S_ref = 195.937,
    a1 = 3.2978e-05,
    a2 = 4796.54,
    a3 = 5.2049e-05,
    a4 = -136106,
    c1 = 109.9262,
    c2 = 186606,
    w_ref = -66110,
    z = 1,
    a_0 = 3.72,
    MM = 0.19836) "Ba(HCO3)+";

  //88
  constant DataRecord Ba_p2(
    R = Modelica.Constants.R/Ba_p2.MM,
    G_ref = -560782,
    H_ref = -537644,
    S_ref = 9.623,
    a1 = 1.1457e-05,
    a2 = -4207.64,
    a3 = -1.966e-06,
    a4 = -98880,
    c1 = 15.8992,
    c2 = -144348,
    w_ref = 412120,
    z = 2,
    a_0 = 5,
    MM = 0.13734) "Ba+2";

  //89
  constant DataRecord BaCl_p1(
    R = Modelica.Constants.R/BaCl_p1.MM,
    G_ref = -689230,
    H_ref = -691711,
    S_ref = 100.416,
    a1 = 1.4449e-05,
    a2 = 273.51,
    a3 = 0.00022954,
    a4 = -117403,
    c1 = 49.4603,
    c2 = -70120,
    w_ref = 79290,
    z = 1,
    a_0 = 3.72,
    MM = 0.17279) "BaCl+";

  //90
  constant DataRecord BaF_p1(
    R = Modelica.Constants.R/BaF_p1.MM,
    G_ref = -841486,
    H_ref = -864038,
    S_ref = 23.012,
    a1 = 4.0384e-06,
    a2 = -2268.48,
    a3 = 0.00032945,
    a4 = -106893,
    c1 = 70.0841,
    c2 = -35685,
    w_ref = 195600,
    z = 1,
    a_0 = 3.72,
    MM = 0.15634) "BaF+";

  //91
  constant DataRecord BaOH_p1(
    R = Modelica.Constants.R/BaOH_p1.MM,
    G_ref = -720903,
    H_ref = -735966,
    S_ref = 115.06,
    a1 = 1.2845e-05,
    a2 = -118.2,
    a3 = 0.00024494,
    a4 = -115780,
    c1 = -23.8756,
    c2 = -317875,
    w_ref = 56990,
    z = 1,
    a_0 = 3.72,
    MM = 0.15435) "BaOH+";

  //92
  constant DataRecord Be_p2(
    R = Modelica.Constants.R/Be_p2.MM,
    G_ref = -378861,
    H_ref = -382836,
    S_ref = -134.306,
    a1 = -2.3924e-06,
    a2 = -3838.49,
    a3 = 0.00039105,
    a4 = -100399,
    c1 = 71.2305,
    c2 = -176397,
    w_ref = 647470,
    z = 2,
    a_0 = 8,
    MM = 0.00901) "Be+2";

  //93
  constant DataRecord BeOH_p1(
    R = Modelica.Constants.R/BeOH_p1.MM,
    G_ref = -585342,
    H_ref = -640989,
    S_ref = -74.894,
    a1 = 9.7496e-06,
    a2 = -876.34,
    a3 = 0.00027527,
    a4 = -112646,
    c1 = 129.5843,
    c2 = 123604,
    w_ref = 344010,
    z = -2,
    a_0 = 3.72,
    MM = 0.026018) "BeOH+";

  //94
  constant DataRecord BeO_aq(
    R = Modelica.Constants.R/BeO_aq.MM,
    G_ref = -538062,
    H_ref = -602914,
    S_ref = -105.018,
    a1 = 7.9149e-06,
    a2 = -1322.52,
    a3 = 0.00029241,
    a4 = -110801,
    c1 = 56.7664,
    c2 = -15313,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.02501) "BeO,aq";

  //95
  constant DataRecord HBeO2_n1(
    R = Modelica.Constants.R/HBeO2_n1.MM,
    G_ref = -720485,
    H_ref = -864833,
    S_ref = -139.327,
    a1 = 1.235e-05,
    a2 = -239.28,
    a3 = 0.00024982,
    a4 = -115282,
    c1 = 235.3655,
    c2 = 315365,
    w_ref = 893330,
    z = -1,
    a_0 = 3.72,
    MM = 0.042018) "HBeO2-";

  //96
  constant DataRecord BeO2_n2(
    R = Modelica.Constants.R/BeO2_n2.MM,
    G_ref = -639734,
    H_ref = -793705,
    S_ref = -171.544,
    a1 = 1.3664e-05,
    a2 = 79.75,
    a3 = 0.0002376,
    a4 = -116600,
    c1 = 142.4091,
    c2 = -235204,
    w_ref = 1603730,
    z = -2,
    a_0 = 3.72,
    MM = 0.04101) "BeO2-2";

  //97
  constant DataRecord Br_n1(
    R = Modelica.Constants.R/Br_n1.MM,
    G_ref = -104056,
    H_ref = -121503,
    S_ref = 82.843,
    a1 = 2.2046e-05,
    a2 = 2758.93,
    a3 = 0.00019853,
    a4 = -131503,
    c1 = -15.8992,
    c2 = -284972,
    w_ref = 579820,
    z = -1,
    a_0 = 3,
    MM = 0.07991) "Br-";

  //98
  constant DataRecord Br3_n1(
    R = Modelica.Constants.R/Br3_n1.MM,
    G_ref = -107069,
    H_ref = -130415,
    S_ref = 215.476,
    a1 = 3.8071e-05,
    a2 = 6040.19,
    a3 = 3.176e-06,
    a4 = -141243,
    c1 = 72.2598,
    c2 = -79241,
    w_ref = 355220,
    z = -1,
    a_0 = 3.72,
    MM = 0.23973) "Br3-";

  //99
  constant DataRecord HBrO_aq(
    R = Modelica.Constants.R/HBrO_aq.MM,
    G_ref = -82425,
    H_ref = -112968,
    S_ref = 141.419,
    a1 = 2.3345e-05,
    a2 = 2443.87,
    a3 = 0.00014462,
    a4 = -126374,
    c1 = 41.2442,
    c2 = -50258,
    w_ref = -71920,
    z = 0,
    a_0 = 3.72,
    MM = 0.096918) "HBrO,aq";

  //100
  constant DataRecord BrO_n1(
    R = Modelica.Constants.R/BrO_n1.MM,
    G_ref = -33472,
    H_ref = -94140,
    S_ref = 40.585,
    a1 = 9.7094e-06,
    a2 = -884.12,
    a3 = 0.00027513,
    a4 = -112612,
    c1 = 14.2691,
    c2 = -365602,
    w_ref = 620070,
    z = -1,
    a_0 = 3.72,
    MM = 0.09591) "BrO-";

  //101
  constant DataRecord BrO3_n1(
    R = Modelica.Constants.R/BrO3_n1.MM,
    G_ref = 18619,
    H_ref = -67070,
    S_ref = 161.712,
    a1 = 2.9128e-05,
    a2 = 3856.52,
    a3 = 8.9002e-05,
    a4 = -132214,
    c1 = 15.5055,
    c2 = -302537,
    w_ref = 436520,
    z = -1,
    a_0 = 3.5,
    MM = 0.12791) "BrO3-";

  //102
  constant DataRecord BrO4_n1(
    R = Modelica.Constants.R/BrO4_n1.MM,
    G_ref = 117989,
    H_ref = 12970,
    S_ref = 199.158,
    a1 = 3.6552e-05,
    a2 = 5669.95,
    a3 = 1.7573e-05,
    a4 = -139708,
    c1 = 48.0675,
    c2 = -171280,
    w_ref = 380120,
    z = -1,
    a_0 = 3.72,
    MM = 0.14391) "BrO4-";

  //103
  constant DataRecord CaOH_p1(
    R = Modelica.Constants.R/CaOH_p1.MM,
    G_ref = -716719,
    H_ref = -751446,
    S_ref = 28.033,
    a1 = 1.1398e-05,
    a2 = -472.92,
    a3 = 0.00025923,
    a4 = -114315,
    c1 = 46.5621,
    c2 = -115031,
    w_ref = 188110,
    z = 1,
    a_0 = 3.72,
    MM = 0.057088) "CaOH+";

  //104
  constant DataRecord Ce_p4(
    R = Modelica.Constants.R/Ce_p4.MM,
    G_ref = -507519,
    H_ref = -576137,
    S_ref = -418.818,
    a1 = -1.7904e-05,
    a2 = -7625.21,
    a3 = 0.00053972,
    a4 = -84747,
    c1 = 168.5273,
    c2 = -126110,
    w_ref = 1546570,
    z = 4,
    a_0 = 11,
    MM = 0) "Ce+4";

  //105
  constant DataRecord HCN_aq(
    R = Modelica.Constants.R/HCN_aq.MM,
    G_ref = 119662,
    H_ref = 107110,
    S_ref = 124.683,
    a1 = 3.3507e-05,
    a2 = 4924.78,
    a3 = 4.7221e-05,
    a4 = -136629,
    c1 = 121.3113,
    c2 = 219911,
    w_ref = -46570,
    z = 0,
    a_0 = 3.72,
    MM = 0.027025) "HCN,aq";

  //106
  constant DataRecord CN_n1(
    R = Modelica.Constants.R/CN_n1.MM,
    G_ref = 172381,
    H_ref = 150624,
    S_ref = 94.14,
    a1 = 2.2892e-05,
    a2 = 2335.22,
    a3 = 0.00014852,
    a4 = -125922,
    c1 = 32.1218,
    c2 = -277818,
    w_ref = 539740,
    z = -1,
    a_0 = 3,
    MM = 0.026017) "CN-";

  //107
  constant DataRecord CO_aq(
    R = Modelica.Constants.R/CO_aq.MM,
    G_ref = -120005,
    H_ref = -120959,
    S_ref = 102.634,
    a1 = 2.6097e-05,
    a2 = 3117,
    a3 = 0.00011792,
    a4 = -129156,
    c1 = 168.4152,
    c2 = 418492,
    w_ref = -155440,
    z = 0,
    a_0 = 3.72,
    MM = 0.02801) "CO,aq";

  //108
  constant DataRecord CO2_aq(
    R = Modelica.Constants.R/CO2_aq.MM,
    G_ref = -385974,
    H_ref = -413798,
    S_ref = 117.57,
    a1 = 2.6136e-05,
    a2 = 3125.91,
    a3 = 0.00011772,
    a4 = -129198,
    c1 = 167.496,
    c2 = 368209,
    w_ref = -8370,
    z = 0,
    a_0 = 3.72,
    MM = 0.04401) "CO2,aq";

  //109
  constant DataRecord CO3_n2(
    R = Modelica.Constants.R/CO3_n2.MM,
    G_ref = -527983,
    H_ref = -675235,
    S_ref = -49.999,
    a1 = 1.1934e-05,
    a2 = -1667.07,
    a3 = 0.00026837,
    a4 = -109382,
    c1 = -13.8934,
    c2 = -719301,
    w_ref = 1418960,
    z = -2,
    a_0 = 4.5,
    MM = 0.06001) "CO3-2";

  //110
  constant DataRecord Ca_CO3_aq(
    R = Modelica.Constants.R/Ca_CO3_aq.MM,
    G_ref = -1099764,
    H_ref = -1202440,
    S_ref = 10.46,
    a1 = -8.786e-07,
    a2 = -3470.63,
    a3 = 0.00037698,
    a4 = -101922,
    c1 = -49.3001,
    c2 = -382920,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.10009) "Ca(CO3),aq";

  //111
  constant DataRecord Ca_HCO3_p1(
    R = Modelica.Constants.R/Ca_HCO3_p1.MM,
    G_ref = -1145705,
    H_ref = -1231560,
    S_ref = 66.944,
    a1 = 1.5506e-05,
    a2 = 530.11,
    a3 = 0.00021974,
    a4 = -118449,
    c1 = 174.5648,
    c2 = 348778,
    w_ref = 128870,
    z = 1,
    a_0 = 3.72,
    MM = 0.1011) "Ca(HCO3)+";

  //112
  constant DataRecord Ca_p2(
    R = Modelica.Constants.R/Ca_p2.MM,
    G_ref = -552790,
    H_ref = -543083,
    S_ref = -56.484,
    a1 = -8.146e-07,
    a2 = -3034.24,
    a3 = 0.00022161,
    a4 = -103730,
    c1 = 37.656,
    c2 = -105520,
    w_ref = 517390,
    z = 2,
    a_0 = 6,
    MM = 0.04008) "Ca+2";

  //113
  constant DataRecord CaCl_p1(
    R = Modelica.Constants.R/CaCl_p1.MM,
    G_ref = -682410,
    H_ref = -705452,
    S_ref = 18.828,
    a1 = 1.1359e-05,
    a2 = -481.03,
    a3 = 0.00025919,
    a4 = -114282,
    c1 = 87.3782,
    c2 = 21928,
    w_ref = 203430,
    z = 1,
    a_0 = 3.72,
    MM = 0.07553) "CaCl+";

  //114
  constant DataRecord CaCl2_aq(
    R = Modelica.Constants.R/CaCl2_aq.MM,
    G_ref = -811696,
    H_ref = -883075,
    S_ref = 25.104,
    a1 = 2.6019e-05,
    a2 = 3098.59,
    a3 = 0.0001185,
    a4 = -129081,
    c1 = 100.2528,
    c2 = 136900,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.11098) "CaCl2,aq";

  //115
  constant DataRecord CaF_p1(
    R = Modelica.Constants.R/CaF_p1.MM,
    G_ref = -838432,
    H_ref = -872782,
    S_ref = -37.656,
    a1 = 6.561e-07,
    a2 = -3094.4,
    a3 = 0.00036191,
    a4 = -103479,
    c1 = 126.2493,
    c2 = 129570,
    w_ref = 289160,
    z = 1,
    a_0 = 3.72,
    MM = 0.05908) "CaF+";

  //116
  constant DataRecord CaSO4_aq(
    R = Modelica.Constants.R/CaSO4_aq.MM,
    G_ref = -1309299,
    H_ref = -1447246,
    S_ref = 20.92,
    a1 = 1.0075e-05,
    a2 = -794.63,
    a3 = 0.00027152,
    a4 = -112985,
    c1 = -35.5397,
    c2 = -340038,
    w_ref = -420,
    z = 0,
    a_0 = 3.72,
    MM = 0.13614) "CaSO4,aq";

  //117
  constant DataRecord Cd_p2(
    R = Modelica.Constants.R/Cd_p2.MM,
    G_ref = -77655,
    H_ref = -75898,
    S_ref = -72.802,
    a1 = 2.247e-07,
    a2 = -4480.23,
    a3 = 0.0006911,
    a4 = -97751,
    c1 = 65.5101,
    c2 = -156800,
    w_ref = 524170,
    z = 2,
    a_0 = 5,
    MM = 0.1124) "Cd+2";

  //118
  constant DataRecord CdOH_p1(
    R = Modelica.Constants.R/CdOH_p1.MM,
    G_ref = -257316,
    H_ref = -307106,
    S_ref = -12.134,
    a1 = 6.163e-07,
    a2 = -3102.06,
    a3 = 0.00036179,
    a4 = -103445,
    c1 = 79.0299,
    c2 = -22133,
    w_ref = 250410,
    z = 1,
    a_0 = 3.72,
    MM = 0.12941) "CdOH+";

  //119
  constant DataRecord CdO_aq(
    R = Modelica.Constants.R/CdO_aq.MM,
    G_ref = -198740,
    H_ref = -249994,
    S_ref = -17.573,
    a1 = -6.155e-07,
    a2 = -3403.31,
    a3 = 0.00037368,
    a4 = -102198,
    c1 = -1.5924,
    c2 = -218158,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1284) "CdO,aq";

  //120
  constant DataRecord HCdO2_n1(
    R = Modelica.Constants.R/HCdO2_n1.MM,
    G_ref = -361916,
    H_ref = -491202,
    S_ref = -45.606,
    a1 = 3.7991e-06,
    a2 = -2325.05,
    a3 = 0.00033127,
    a4 = -106659,
    c1 = 97.2554,
    c2 = -119294,
    w_ref = 751660,
    z = -1,
    a_0 = 3.72,
    MM = 0.14541) "HCdO2-";

  //121
  constant DataRecord CdO2_n2(
    R = Modelica.Constants.R/CdO2_n2.MM,
    G_ref = -281583,
    H_ref = -422584,
    S_ref = -84.935,
    a1 = 3.3146e-06,
    a2 = -2443.5,
    a3 = 0.00033599,
    a4 = -106169,
    c1 = 14.5645,
    c2 = -637478,
    w_ref = 1472350,
    z = -2,
    a_0 = 3.72,
    MM = 0.1444) "CdO2-2";

  //122
  constant DataRecord Ce_p2(
    R = Modelica.Constants.R/Ce_p2.MM,
    G_ref = -313382,
    H_ref = -294135,
    S_ref = 5.858,
    a1 = 5.502e-07,
    a2 = -3118.67,
    a3 = 0.00036256,
    a4 = -103378,
    c1 = 35.8799,
    c2 = -230940,
    w_ref = 434130,
    z = 2,
    a_0 = 3.72,
    MM = 0.14012) "Ce+2";

  //123
  constant DataRecord Ce_p3(
    R = Modelica.Constants.R/Ce_p3.MM,
    G_ref = -676134,
    H_ref = -700402,
    S_ref = -205.016,
    a1 = -1.4574e-05,
    a2 = -6811.09,
    a3 = 0.00050753,
    a4 = -88111,
    c1 = -1.4853,
    c2 = -533502,
    w_ref = 973410,
    z = 3,
    a_0 = 9,
    MM = 0.14012) "Ce+3";

  //124
  constant DataRecord Cl_n1(
    R = Modelica.Constants.R/Cl_n1.MM,
    G_ref = -131290,
    H_ref = -167080,
    S_ref = 56.735,
    a1 = 1.687e-05,
    a2 = 2008.74,
    a3 = 0.00023276,
    a4 = -119118,
    c1 = -18.4096,
    c2 = -239074,
    w_ref = 609190,
    z = -1,
    a_0 = 3,
    MM = 0.03545) "Cl-";

  //125
  constant DataRecord HClO_aq(
    R = Modelica.Constants.R/HClO_aq.MM,
    G_ref = -79914,
    H_ref = -120918,
    S_ref = 141.838,
    a1 = 2.34e-05,
    a2 = 2458.14,
    a3 = 0.00014388,
    a4 = -126432,
    c1 = 41.4312,
    c2 = -49405,
    w_ref = -72550,
    z = 0,
    a_0 = 3.72,
    MM = 0.052458) "HClO,aq";

  //126
  constant DataRecord ClO_n1(
    R = Modelica.Constants.R/ClO_n1.MM,
    G_ref = -36819,
    H_ref = -107110,
    S_ref = 41.84,
    a1 = 9.8738e-06,
    a2 = -843.66,
    a3 = 0.00027345,
    a4 = -112780,
    c1 = 14.5549,
    c2 = -363899,
    w_ref = 617850,
    z = -1,
    a_0 = 3.72,
    MM = 0.05145) "ClO-";

  //127
  constant DataRecord ClO2_n1(
    R = Modelica.Constants.R/ClO2_n1.MM,
    G_ref = 17154,
    H_ref = -66526,
    S_ref = 101.253,
    a1 = 2.1825e-05,
    a2 = 2074.43,
    a3 = 0.00015878,
    a4 = -124846,
    c1 = 27.1864,
    c2 = -291453,
    w_ref = 528730,
    z = -1,
    a_0 = 4.25,
    MM = 0.06745) "ClO2-";

  //128
  constant DataRecord ClO3_n1(
    R = Modelica.Constants.R/ClO3_n1.MM,
    G_ref = -7950,
    H_ref = -103972,
    S_ref = 162.339,
    a1 = 2.9985e-05,
    a2 = 4065.68,
    a3 = 8.078e-05,
    a4 = -133080,
    c1 = 35.7987,
    c2 = -231798,
    w_ref = 435890,
    z = -1,
    a_0 = 3.5,
    MM = 0.08345) "ClO3-";

  //129
  constant DataRecord ClO4_n1(
    R = Modelica.Constants.R/ClO4_n1.MM,
    G_ref = -8535,
    H_ref = -129327,
    S_ref = 182.004,
    a1 = 3.4062e-05,
    a2 = 6512.56,
    a3 = -0.00032667,
    a4 = -143193,
    c1 = 68.8268,
    c2 = -274889,
    w_ref = 405810,
    z = -1,
    a_0 = 3.5,
    MM = 0.09945) "ClO4-";

  //130
  constant DataRecord Co_p2(
    R = Modelica.Constants.R/Co_p2.MM,
    G_ref = -54392,
    H_ref = -58158,
    S_ref = -112.968,
    a1 = -5.1262e-06,
    a2 = -3738.66,
    a3 = 0.00022255,
    a4 = -100813,
    c1 = 63.6022,
    c2 = -193443,
    w_ref = 617930,
    z = 2,
    a_0 = 6,
    MM = 0.05893) "Co+2";

  //131
  constant DataRecord Co_p3(
    R = Modelica.Constants.R/Co_p3.MM,
    G_ref = 133888,
    H_ref = 92048,
    S_ref = -305.432,
    a1 = -1.1999e-05,
    a2 = -6182.99,
    a3 = 0.000483,
    a4 = -90709,
    c1 = 68.9268,
    c2 = -337477,
    w_ref = 1125540,
    z = 3,
    a_0 = 3.72,
    MM = 0.05893) "Co+3";

  //132
  constant DataRecord CoOH_p1(
    R = Modelica.Constants.R/CoOH_p1.MM,
    G_ref = -234409,
    H_ref = -286395,
    S_ref = -41.84,
    a1 = -1.0715e-06,
    a2 = -3516.19,
    a3 = 0.00037847,
    a4 = -101734,
    c1 = 111.6454,
    c2 = 77580,
    w_ref = 293010,
    z = 1,
    a_0 = 3.72,
    MM = 0.075938) "CoOH+";

  //133
  constant DataRecord CoO_aq(
    R = Modelica.Constants.R/CoO_aq.MM,
    G_ref = -184096,
    H_ref = -246019,
    S_ref = -74.475,
    a1 = -6.5128e-06,
    a2 = -4842.39,
    a3 = 0.00043008,
    a4 = -96249,
    c1 = 36.4138,
    c2 = -86052,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.07493) "CoO,aq";

  //134
  constant DataRecord HCoO2_n1(
    R = Modelica.Constants.R/HCoO2_n1.MM,
    G_ref = -348946,
    H_ref = -489946,
    S_ref = -106.692,
    a1 = -4.063e-07,
    a2 = -3353.69,
    a3 = 0.00037207,
    a4 = -102403,
    c1 = 187.3227,
    c2 = 163661,
    w_ref = 845630,
    z = -1,
    a_0 = 3.72,
    MM = 0.091938) "HCoO2-";

  //135
  constant DataRecord CoO2_n2(
    R = Modelica.Constants.R/CoO2_n2.MM,
    G_ref = -264429,
    H_ref = -414216,
    S_ref = -136.817,
    a1 = -2.469e-07,
    a2 = -3314.52,
    a3 = 0.00037045,
    a4 = -102567,
    c1 = 91.4279,
    c2 = -396283,
    w_ref = 1553440,
    z = -2,
    a_0 = 3.72,
    MM = 0.09093) "CoO2-2";

  //136
  constant DataRecord CoOH_p2(
    R = Modelica.Constants.R/CoOH_p2.MM,
    G_ref = -96232,
    H_ref = -151042,
    S_ref = -116.315,
    a1 = -4.5539e-06,
    a2 = -4366.13,
    a3 = 0.00041181,
    a4 = -98219,
    c1 = 69.2419,
    c2 = -173841,
    w_ref = 617930,
    z = 2,
    a_0 = 3.72,
    MM = 0.075938) "CoOH+2";

  //137
  constant DataRecord o_CRESOL_aq(
    R = Modelica.Constants.R/o_CRESOL_aq.MM,
    G_ref = -53241,
    H_ref = -188619,
    S_ref = 210.874,
    a1 = 6.4509e-05,
    a2 = 11806.58,
    a3 = -7.5643e-05,
    a4 = -165076,
    c1 = 376.1788,
    c2 = 217384,
    w_ref = -105940,
    z = 0,
    a_0 = 3.72,
    MM = 0.10813) "o-CRESOL,aq";

  //138
  constant DataRecord m_CRESOL_aq(
    R = Modelica.Constants.R/m_CRESOL_aq.MM,
    G_ref = -58116,
    H_ref = -190623,
    S_ref = 220.497,
    a1 = 6.4531e-05,
    a2 = 11811.35,
    a3 = -7.5735e-05,
    a4 = -165096,
    c1 = 371.7505,
    c2 = 211183,
    w_ref = -99580,
    z = 0,
    a_0 = 3.72,
    MM = 0.10813) "m-CRESOL,aq";

  //139
  constant DataRecord p_CRESOL_aq(
    R = Modelica.Constants.R/p_CRESOL_aq.MM,
    G_ref = -48639,
    H_ref = -185121,
    S_ref = 207.108,
    a1 = 6.4501e-05,
    a2 = 11806.16,
    a3 = -7.5977e-05,
    a4 = -165076,
    c1 = 339.1651,
    c2 = 171912,
    w_ref = -108450,
    z = 0,
    a_0 = 3.72,
    MM = 0.10813) "p-CRESOL,aq";

  //140
  constant DataRecord Cr2O7_n2(
    R = Modelica.Constants.R/Cr2O7_n2.MM,
    G_ref = -1309466,
    H_ref = -1486868,
    S_ref = 301.248,
    a1 = 5.2008e-05,
    a2 = 9442.45,
    a3 = -0.00013038,
    a4 = -155306,
    c1 = 36.2251,
    c2 = -374978,
    w_ref = 887680,
    z = -2,
    a_0 = 3.72,
    MM = 0.216) "Cr2O7-2";

  //141
  constant DataRecord CrO4_n2(
    R = Modelica.Constants.R/CrO4_n2.MM,
    G_ref = -731363,
    H_ref = -882531,
    S_ref = 57.739,
    a1 = 2.2966e-05,
    a2 = 2352.37,
    a3 = 0.00014804,
    a4 = -125993,
    c1 = -11.9696,
    c2 = -660490,
    w_ref = 1256200,
    z = -2,
    a_0 = 4,
    MM = 0.116) "CrO4-2";

  //142
  constant DataRecord Cs_p1(
    R = Modelica.Constants.R/Cs_p1.MM,
    G_ref = -291667,
    H_ref = -258027,
    S_ref = 132.842,
    a1 = 2.5721e-05,
    a2 = -54.77,
    a3 = 0.00017612,
    a4 = -116047,
    c1 = 26.2337,
    c2 = -239994,
    w_ref = 40750,
    z = 1,
    a_0 = 3.72,
    MM = 0.13291) "Cs+";

  //143
  constant DataRecord CsBr_aq(
    R = Modelica.Constants.R/CsBr_aq.MM,
    G_ref = -395848,
    H_ref = -370242,
    S_ref = 239.743,
    a1 = 4.0018e-05,
    a2 = 6516.87,
    a3 = -1.5853e-05,
    a4 = -143210,
    c1 = -39.3894,
    c2 = -353418,
    w_ref = -420,
    z = 0,
    a_0 = 3.72,
    MM = 0.21282) "CsBr,aq";

  //144
  constant DataRecord CsCl_aq(
    R = Modelica.Constants.R/CsCl_aq.MM,
    G_ref = -422166,
    H_ref = -422375,
    S_ref = 219.995,
    a1 = 3.515e-05,
    a2 = 5328.07,
    a3 = 3.0874e-05,
    a4 = -138298,
    c1 = -25.3396,
    c2 = -331515,
    w_ref = 83680,
    z = 0,
    a_0 = 3.72,
    MM = 0.16836) "CsCl,aq";

  //145
  constant DataRecord CsI_aq(
    R = Modelica.Constants.R/CsI_aq.MM,
    G_ref = -349197,
    H_ref = -325599,
    S_ref = 252.295,
    a1 = 4.8749e-05,
    a2 = 8648.7,
    a3 = -9.9642e-05,
    a4 = -152026,
    c1 = -30.175,
    c2 = -334921,
    w_ref = 41840,
    z = 0,
    a_0 = 3.72,
    MM = 0.25981) "CsI,aq";

  //146
  constant DataRecord CsOH_aq(
    R = Modelica.Constants.R/CsOH_aq.MM,
    G_ref = -439320,
    H_ref = -469863,
    S_ref = 150.206,
    a1 = 2.3717e-05,
    a2 = 2534.29,
    a3 = 0.00014117,
    a4 = -126746,
    c1 = -53.577,
    c2 = -398840,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14992) "CsOH,aq";

  //147
  constant DataRecord Cu_p1(
    R = Modelica.Constants.R/Cu_p1.MM,
    G_ref = 49999,
    H_ref = 71680,
    S_ref = 40.585,
    a1 = 3.3765e-06,
    a2 = -2428.39,
    a3 = 0.00033541,
    a4 = -106232,
    c1 = 74.9911,
    c2 = -10201,
    w_ref = 169280,
    z = 1,
    a_0 = 3.72,
    MM = 0.06354) "Cu+";

  //148
  constant DataRecord Cu_p2(
    R = Modelica.Constants.R/Cu_p2.MM,
    G_ref = 65584,
    H_ref = 65689,
    S_ref = -97.069,
    a1 = -4.6112e-06,
    a2 = -4381.74,
    a3 = 0.0004128,
    a4 = -98161,
    c1 = 84.9352,
    c2 = -183678,
    w_ref = 617930,
    z = 2,
    a_0 = 6,
    MM = 0.06354) "Cu+2";

  //149
  constant DataRecord CuOH_p1(
    R = Modelica.Constants.R/CuOH_p1.MM,
    G_ref = -126357,
    H_ref = -174473,
    S_ref = -25.522,
    a1 = -1.159e-07,
    a2 = -3280.93,
    a3 = 0.00036879,
    a4 = -102705,
    c1 = 89.9819,
    c2 = 9397,
    w_ref = 270790,
    z = 1,
    a_0 = 3.72,
    MM = 0.080548) "CuOH+";

  //150
  constant DataRecord CuO_aq(
    R = Modelica.Constants.R/CuO_aq.MM,
    G_ref = -87027,
    H_ref = -143093,
    S_ref = -51.882,
    a1 = -1.6472e-06,
    a2 = -3655.14,
    a3 = 0.00038357,
    a4 = -101161,
    c1 = 50.8314,
    c2 = -138896,
    w_ref = 308950,
    z = 0,
    a_0 = 3.72,
    MM = 0.07954) "CuO,aq";

  //151
  constant DataRecord HCuO2_n1(
    R = Modelica.Constants.R/HCuO2_n1.MM,
    G_ref = -251458,
    H_ref = -361079,
    S_ref = 1.255,
    a1 = 5.5237e-06,
    a2 = -1904.85,
    a3 = 0.00031498,
    a4 = -108395,
    c1 = -13.7022,
    c2 = -337477,
    w_ref = 228610,
    z = -1,
    a_0 = 3.72,
    MM = 0.096548) "HCuO2-";

  //152
  constant DataRecord CuO2_n2(
    R = Modelica.Constants.R/CuO2_n2.MM,
    G_ref = -172381,
    H_ref = -311290,
    S_ref = -96.65,
    a1 = 2.5204e-06,
    a2 = -2637.51,
    a3 = 0.00034363,
    a4 = -105366,
    c1 = 32.021,
    c2 = -582936,
    w_ref = 1491470,
    z = -2,
    a_0 = 3.72,
    MM = 0.09554) "CuO2-2";

  //153
  constant DataRecord DECANOATE_aq(
    R = Modelica.Constants.R/DECANOATE_aq.MM,
    G_ref = -303173,
    H_ref = -680737,
    S_ref = 301.666,
    a1 = 0.00010272,
    a2 = 20355.75,
    a3 = -0.00024407,
    a4 = -200418,
    c1 = 759.9818,
    c2 = -75730,
    w_ref = 224760,
    z = -1,
    a_0 = 3.72,
    MM = 0.17125) "DECANOATE,aq";

  //154
  constant DataRecord DECANOIC_ACID_aq(
    R = Modelica.Constants.R/DECANOIC_ACID_aq.MM,
    G_ref = -331247,
    H_ref = -678645,
    S_ref = 402.919,
    a1 = 0.00010994,
    a2 = 22103.82,
    a3 = -0.00030962,
    a4 = -207648,
    c1 = 745.7135,
    c2 = 659700,
    w_ref = 21340,
    z = 0,
    a_0 = 3.72,
    MM = 0.17226) "DECANOIC-ACID,aq";

  //155
  constant DataRecord TWO_3DMP_aq(
    R = Modelica.Constants.R/TWO_3DMP_aq.MM,
    G_ref = -46045,
    H_ref = -225304,
    S_ref = 199.995,
    a1 = 7.2443e-05,
    a2 = 13605.49,
    a3 = -0.00011668,
    a4 = -172515,
    c1 = 443.3965,
    c2 = 301302,
    w_ref = -113180,
    z = 0,
    a_0 = 3.72,
    MM = 0.12216) "2-3DMP,aq";

  //156
  constant DataRecord TWO_4DMP_aq(
    R = Modelica.Constants.R/TWO_4DMP_aq.MM,
    G_ref = -54731,
    H_ref = -220392,
    S_ref = 245.601,
    a1 = 7.2545e-05,
    a2 = 13628.33,
    a3 = -0.00011713,
    a4 = -172607,
    c1 = 428.4596,
    c2 = 279391,
    w_ref = -82930,
    z = 0,
    a_0 = 3.72,
    MM = 0.12216) "2,4-DMP,aq";

  //157
  constant DataRecord TWO_5DMP_aq(
    R = Modelica.Constants.R/TWO_5DMP_aq.MM,
    G_ref = -49685,
    H_ref = -227693,
    S_ref = 204.179,
    a1 = 7.2453e-05,
    a2 = 13608.29,
    a3 = -0.00011688,
    a4 = -172527,
    c1 = 463.7156,
    c2 = 326105,
    w_ref = -110420,
    z = 0,
    a_0 = 3.72,
    MM = 0.12216) "2-5DMP,aq";

  //158
  constant DataRecord TWO_6DMP_aq(
    R = Modelica.Constants.R/TWO_6DMP_aq.MM,
    G_ref = -45911,
    H_ref = -219430,
    S_ref = 219.242,
    a1 = 7.2487e-05,
    a2 = 13615.07,
    a3 = -0.00011686,
    a4 = -172552,
    c1 = 430.1934,
    c2 = 283525,
    w_ref = -100420,
    z = 0,
    a_0 = 3.72,
    MM = 0.12216) "2-6DMP,aq";

  //159
  constant DataRecord THREE_4DMP_aq(
    R = Modelica.Constants.R/THREE_4DMP_aq.MM,
    G_ref = -51313,
    H_ref = -227702,
    S_ref = 209.618,
    a1 = 7.2465e-05,
    a2 = 13610.22,
    a3 = -0.00011676,
    a4 = -172536,
    c1 = 409.8759,
    c2 = 259136,
    w_ref = -106780,
    z = 0,
    a_0 = 3.72,
    MM = 0.12216) "3-4DMP,aq";

  //160
  constant DataRecord THREE_5DMP_aq(
    R = Modelica.Constants.R/THREE_5DMP_aq.MM,
    G_ref = -52128,
    H_ref = -229890,
    S_ref = 205.016,
    a1 = 7.2455e-05,
    a2 = 13607.75,
    a3 = -0.00011667,
    a4 = -172523,
    c1 = 437.6845,
    c2 = 293859,
    w_ref = -109830,
    z = 0,
    a_0 = 3.72,
    MM = 0.12216) "3-5DMP,aq";

  //161
  constant DataRecord DODECANOATE_aq(
    R = Modelica.Constants.R/DODECANOATE_aq.MM,
    G_ref = -286060,
    H_ref = 728183,
    S_ref = 357.732,
    a1 = 0.00012013,
    a2 = 24413.43,
    a3 = -0.00036209,
    a4 = -217196,
    c1 = 908.0481,
    c2 = -64894,
    w_ref = 139750,
    z = -1,
    a_0 = 3.72,
    MM = 0.1993) "DODECANOATE,aq";

  //162
  constant DataRecord DODECANOIC_ACID_aq(
    R = Modelica.Constants.R/DODECANOIC_ACID_aq.MM,
    G_ref = -314135,
    H_ref = -726091,
    S_ref = 458.985,
    a1 = 0.00012816,
    a2 = 26232.26,
    a3 = -0.00040328,
    a4 = -224714,
    c1 = 890.9217,
    c2 = 834976,
    w_ref = 58530,
    z = 0,
    a_0 = 3.72,
    MM = 0.20031) "DODECANOIC-ACID,aq";

  //163
  constant DataRecord Dy_p2(
    R = Modelica.Constants.R/Dy_p2.MM,
    G_ref = -430115,
    H_ref = -417982,
    S_ref = -15.062,
    a1 = 8.62e-08,
    a2 = -3232.06,
    a3 = 0.00036704,
    a4 = -102910,
    c1 = 41.2919,
    c2 = -222417,
    w_ref = 466260,
    z = 2,
    a_0 = 3.72,
    MM = 0.1625) "Dy+2";

  //164
  constant DataRecord Dy_p3(
    R = Modelica.Constants.R/Dy_p3.MM,
    G_ref = -664001,
    H_ref = -696636,
    S_ref = -230.957,
    a1 = -1.2553e-05,
    a2 = -6320.94,
    a3 = 0.00048902,
    a4 = -90144,
    c1 = 39.7798,
    c2 = -397141,
    w_ref = 995460,
    z = 3,
    a_0 = 3.72,
    MM = 0.1625) "Dy+3";

  //165
  constant DataRecord Dy_p4(
    R = Modelica.Constants.R/Dy_p4.MM,
    G_ref = -233049,
    H_ref = -307106,
    S_ref = -435.136,
    a1 = -1.806e-05,
    a2 = -7662.74,
    a3 = 0.00054112,
    a4 = -84592,
    c1 = 167.5901,
    c2 = -136340,
    w_ref = 1568330,
    z = 4,
    a_0 = 3.72,
    MM = 0.1625) "Dy+4";

  //166
  constant DataRecord ETHANAMINE_aq(
    R = Modelica.Constants.R/ETHANAMINE_aq.MM,
    G_ref = 26359,
    H_ref = -99705,
    S_ref = 140.582,
    a1 = 4.0702e-05,
    a2 = 4738.67,
    a3 = 0.00047133,
    a4 = -135859,
    c1 = 232.9647,
    c2 = -8962,
    w_ref = -70630,
    z = 0,
    a_0 = 3.72,
    MM = 0.045082) "ETHANAMINE,aq";

  //167
  constant DataRecord ETHANE_aq(
    R = Modelica.Constants.R/ETHANE_aq.MM,
    G_ref = -16259,
    H_ref = -103387,
    S_ref = 111.253,
    a1 = 3.661e-05,
    a2 = 5483.17,
    a3 = 6.8023e-05,
    a4 = -138938,
    c1 = 308.8093,
    c2 = 93061,
    w_ref = -26230,
    z = 0,
    a_0 = 3.72,
    MM = 0.030067) "ETHANE,aq";

  //168
  constant DataRecord ETHANOL_aq(
    R = Modelica.Constants.R/ETHANOL_aq.MM,
    G_ref = -181293,
    H_ref = -287232,
    S_ref = 150.206,
    a1 = 3.8632e-05,
    a2 = 4166.47,
    a3 = 0.00050813,
    a4 = -133495,
    c1 = 251.1132,
    c2 = 6305,
    w_ref = -85230,
    z = 0,
    a_0 = 3.72,
    MM = 0.046067) "ETHANOL,aq";

  //169
  constant DataRecord ETHYLACETATE_aq(
    R = Modelica.Constants.R/ETHYLACETATE_aq.MM,
    G_ref = -331791,
    H_ref = -488859,
    S_ref = 223.426,
    a1 = 5.7894e-06,
    a2 = 10309.21,
    a3 = -4.1978e-05,
    a4 = -158887,
    c1 = 350.1937,
    c2 = 184314,
    w_ref = -97650,
    z = 0,
    a_0 = 3.72,
    MM = 0.088103) "ETHYLACETATE,aq";

  //170
  constant DataRecord ETHYLBENZENE_aq(
    R = Modelica.Constants.R/ETHYLBENZENE_aq.MM,
    G_ref = 135729,
    H_ref = -10460,
    S_ref = 208.363,
    a1 = 7.3459e-05,
    a2 = 13835.44,
    a3 = -0.00012186,
    a4 = -173464,
    c1 = 429.06,
    c2 = 280219,
    w_ref = -83680,
    z = 0,
    a_0 = 3.72,
    MM = 0.10616) "ETHYLBENZENE,aq";

  //171
  constant DataRecord ETHYLENE_aq(
    R = Modelica.Constants.R/ETHYLENE_aq.MM,
    G_ref = 81379,
    H_ref = 35857,
    S_ref = 120.081,
    a1 = 3.287e-05,
    a2 = 5288.2,
    a3 = -7.8396e-05,
    a4 = -138131,
    c1 = 163.5944,
    c2 = 405848,
    w_ref = -167360,
    z = 0,
    a_0 = 3.72,
    MM = 0.028052) "ETHYLENE,aq";

  //172
  constant DataRecord ETHYNE_aq(
    R = Modelica.Constants.R/ETHYNE_aq.MM,
    G_ref = 217108,
    H_ref = 212129,
    S_ref = 125.52,
    a1 = 3.4969e-05,
    a2 = 5283.22,
    a3 = 3.2773e-05,
    a4 = -138110,
    c1 = 138.249,
    c2 = 324741,
    w_ref = -190080,
    z = 0,
    a_0 = 3.72,
    MM = 0.026036) "ETHYNE,aq";

  //173
  constant DataRecord Er_p2(
    R = Modelica.Constants.R/Er_p2.MM,
    G_ref = -383254,
    H_ref = -373004,
    S_ref = -23.012,
    a1 = -4.35e-08,
    a2 = -3263.35,
    a3 = 0.00036811,
    a4 = -102780,
    c1 = 43.158,
    c2 = -219861,
    w_ref = 478520,
    z = 2,
    a_0 = 3.72,
    MM = 0.16726) "Er+2";

  //174
  constant DataRecord Er_p3(
    R = Modelica.Constants.R/Er_p3.MM,
    G_ref = -669022,
    H_ref = -705004,
    S_ref = -243.927,
    a1 = -1.3824e-05,
    a2 = -6631.31,
    a3 = 0.00050122,
    a4 = -88860,
    c1 = 34.6498,
    c2 = -419300,
    w_ref = 1008970,
    z = 3,
    a_0 = 3.72,
    MM = 0.16726) "Er+3";

  //175
  constant DataRecord Er_p4(
    R = Modelica.Constants.R/Er_p4.MM,
    G_ref = -117570,
    H_ref = -192046,
    S_ref = -438.065,
    a1 = -1.8098e-05,
    a2 = -7673.54,
    a3 = 0.00054181,
    a4 = -84546,
    c1 = 167.3613,
    c2 = -138896,
    w_ref = 1573810,
    z = 4,
    a_0 = 3.72,
    MM = 0.16726) "Er+4";

  //176
  constant DataRecord Eu_p2(
    R = Modelica.Constants.R/Eu_p2.MM,
    G_ref = -540573,
    H_ref = -527812,
    S_ref = -10.042,
    a1 = 1.703e-07,
    a2 = -3212.31,
    a3 = 0.00036643,
    a4 = -102989,
    c1 = 39.9735,
    c2 = -224124,
    w_ref = 457270,
    z = 2,
    a_0 = 3.72,
    MM = 0.15196) "Eu+2";

  //177
  constant DataRecord Eu_p3(
    R = Modelica.Constants.R/Eu_p3.MM,
    G_ref = -574463,
    H_ref = -605425,
    S_ref = -221.752,
    a1 = -1.2986e-05,
    a2 = -6426.58,
    a3 = 0.00049317,
    a4 = -89705,
    c1 = 25.3333,
    c2 = -438902,
    w_ref = 969060,
    z = 3,
    a_0 = 3.72,
    MM = 0.15196) "Eu+3";

  //178
  constant DataRecord Eu_p4(
    R = Modelica.Constants.R/Eu_p4.MM,
    G_ref = 20920,
    H_ref = -53137,
    S_ref = -432.207,
    a1 = -1.8021e-05,
    a2 = -7654.67,
    a3 = 0.00054108,
    a4 = -84626,
    c1 = 167.5755,
    c2 = -134633,
    w_ref = 1562850,
    z = 4,
    a_0 = 3.72,
    MM = 0.15196) "Eu+4";

  //179
  constant DataRecord F_n1(
    R = Modelica.Constants.R/F_n1.MM,
    G_ref = -281751,
    H_ref = -335348,
    S_ref = -13.18,
    a1 = 2.8744e-06,
    a2 = 568.52,
    a3 = 0.00031812,
    a4 = -118625,
    c1 = 18.6606,
    c2 = -313298,
    w_ref = 747680,
    z = -1,
    a_0 = 3.5,
    MM = 0.019) "F-";

  //180
  constant DataRecord FORMATE_aq(
    R = Modelica.Constants.R/FORMATE_aq.MM,
    G_ref = -350879,
    H_ref = -425429,
    S_ref = 90.793,
    a1 = 2.4201e-05,
    a2 = 1976.61,
    a3 = 0.00030807,
    a4 = -124441,
    c1 = 71.128,
    c2 = -518816,
    w_ref = 544050,
    z = -1,
    a_0 = 3.72,
    MM = 0.045018) "FORMATE,aq";

  //181
  constant DataRecord FORMIC_ACID_aq(
    R = Modelica.Constants.R/FORMIC_ACID_aq.MM,
    G_ref = -372301,
    H_ref = -425429,
    S_ref = 162.758,
    a1 = 2.676e-05,
    a2 = 3251.51,
    a3 = 0.00011848,
    a4 = -129712,
    c1 = 109.2024,
    c2 = -129704,
    w_ref = -138070,
    z = 0,
    a_0 = 3.72,
    MM = 0.046026) "FORMIC-ACID,aq";

  //182
  constant DataRecord Fe_p2(
    R = Modelica.Constants.R/Fe_p2.MM,
    G_ref = -91504,
    H_ref = -92257,
    S_ref = -105.855,
    a1 = -3.2916e-06,
    a2 = -4057.18,
    a3 = 0.00039948,
    a4 = -99496,
    c1 = 61.8646,
    c2 = -194292,
    w_ref = 601740,
    z = 2,
    a_0 = 6,
    MM = 0.05585) "Fe+2";

  //183
  constant DataRecord Fe_p3(
    R = Modelica.Constants.R/Fe_p3.MM,
    G_ref = -17238,
    H_ref = -49580,
    S_ref = -277.399,
    a1 = -1.0149e-05,
    a2 = -5730.45,
    a3 = 0.00046501,
    a4 = -92579,
    c1 = 79.688,
    c2 = -285487,
    w_ref = 1079970,
    z = 3,
    a_0 = 9,
    MM = 0.05585) "Fe+3";

  //184
  constant DataRecord FeCl_p1(
    R = Modelica.Constants.R/FeCl_p1.MM,
    G_ref = -221878,
    H_ref = -256312,
    S_ref = -42.091,
    a1 = 8.9822e-06,
    a2 = -1061.36,
    a3 = 0.00028201,
    a4 = -111884,
    c1 = 103.308,
    c2 = 48606,
    w_ref = 293010,
    z = 1,
    a_0 = 3.72,
    MM = 0.0913) "FeCl+";

  //185
  constant DataRecord FeCl2_aq(
    R = Modelica.Constants.R/FeCl2_aq.MM,
    G_ref = -307440,
    H_ref = -328402,
    S_ref = 179.912,
    a1 = 2.3036e-05,
    a2 = 2370.24,
    a3 = 0.00014713,
    a4 = -126068,
    c1 = 95.937,
    c2 = 121901,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.12675) "FeCl2,aq";

  //186
  constant DataRecord FeOH_p2(
    R = Modelica.Constants.R/FeOH_p2.MM,
    G_ref = -241835,
    H_ref = -292880,
    S_ref = -106.274,
    a1 = -4.8375e-06,
    a2 = -4435.42,
    a3 = 0.00041454,
    a4 = -97935,
    c1 = 61.1291,
    c2 = -196849,
    w_ref = 601740,
    z = 2,
    a_0 = 3.72,
    MM = 0.072858) "FeOH+2";

  //187
  constant DataRecord FeOH_p1(
    R = Modelica.Constants.R/FeOH_p1.MM,
    G_ref = -275516,
    H_ref = -326687,
    S_ref = -41.84,
    a1 = -1.0715e-06,
    a2 = -3516.19,
    a3 = 0.00037847,
    a4 = -101734,
    c1 = 89.5765,
    c2 = 874,
    w_ref = 293010,
    z = 1,
    a_0 = 3.72,
    MM = 0.072858) "FeOH+";

  //188
  constant DataRecord FeO_p1(
    R = Modelica.Constants.R/FeO_p1.MM,
    G_ref = -222170,
    H_ref = -255224,
    S_ref = -46.442,
    a1 = -1.553e-05,
    a2 = -7046.19,
    a3 = 0.00051712,
    a4 = -87140,
    c1 = -64.4261,
    c2 = -536912,
    w_ref = 300870,
    z = 1,
    a_0 = 3.72,
    MM = 0.07185) "FeO+";

  //189
  constant DataRecord FeO(
    R = Modelica.Constants.R/FeO.MM,
    G_ref = -212212,
    H_ref = -263383,
    S_ref = -41.84,
    a1 = -2.1041e-06,
    a2 = -3767.82,
    a3 = 0.00038824,
    a4 = -100692,
    c1 = 24.6442,
    c2 = -126963,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.07185) "FeO";

  //190
  constant DataRecord HFeO2_n1(
    R = Modelica.Constants.R/HFeO2_n1.MM,
    G_ref = -399154,
    H_ref = -525929,
    S_ref = -62.76,
    a1 = 2.6242e-06,
    a2 = -2612.49,
    a3 = 0.00034269,
    a4 = -105470,
    c1 = 151.3035,
    c2 = 60534,
    w_ref = 776720,
    z = -1,
    a_0 = 3.72,
    MM = 0.088858) "HFeO2-";

  //191
  constant DataRecord HFeO2_aq(
    R = Modelica.Constants.R/HFeO2_aq.MM,
    G_ref = -423002,
    H_ref = -503335,
    S_ref = 92.885,
    a1 = 1.1465e-05,
    a2 = -456.27,
    a3 = 0.00025847,
    a4 = -114382,
    c1 = -158.2807,
    c2 = -762764,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.088858) "HFeO2,aq";

  //192
  constant DataRecord FeO2_n1(
    R = Modelica.Constants.R/FeO2_n1.MM,
    G_ref = -368192,
    H_ref = -463169,
    S_ref = 44.35,
    a1 = 9.9734e-06,
    a2 = -820.15,
    a3 = 0.00027272,
    a4 = -112880,
    c1 = -55.7338,
    c2 = -606797,
    w_ref = 613460,
    z = -1,
    a_0 = 3.72,
    MM = 0.08785) "FeO2-";

  //193
  constant DataRecord GLUTAMIC_ACID_aq(
    R = Modelica.Constants.R/GLUTAMIC_ACID_aq.MM,
    G_ref = -723832,
    H_ref = -970688,
    S_ref = 294.972,
    a1 = 5.7518e-05,
    a2 = 8100.48,
    a3 = 0.0004988,
    a4 = -149758,
    c1 = 159.0903,
    c2 = -49974,
    w_ref = -304430,
    z = 0,
    a_0 = 3.72,
    MM = 0.14713) "GLUTAMIC-ACID,aq";

  //194
  constant DataRecord GLUTAMINE_aq(
    R = Modelica.Constants.R/GLUTAMINE_aq.MM,
    G_ref = -529694,
    H_ref = -804709,
    S_ref = 258.99,
    a1 = 6.0565e-05,
    a2 = 8709.96,
    a3 = 0.0005037,
    a4 = -152277,
    c1 = 172.7733,
    c2 = -43187,
    w_ref = -249950,
    z = 0,
    a_0 = 3.72,
    MM = 0.14614) "GLUTAMINE,aq";

  //195
  constant DataRecord GLUTARATE_aq(
    R = Modelica.Constants.R/GLUTARATE_aq.MM,
    G_ref = -683958,
    H_ref = -937802,
    S_ref = 109.621,
    a1 = 5.2704e-05,
    a2 = 9131.5,
    a3 = -1.5029e-05,
    a4 = -154017,
    c1 = -18.5602,
    c2 = -191924,
    w_ref = 1178170,
    z = -2,
    a_0 = 3.72,
    MM = 0.1301) "GLUTARATE,aq";

  //196
  constant DataRecord GLUTARIC_ACID_aq(
    R = Modelica.Constants.R/GLUTARIC_ACID_aq.MM,
    G_ref = -739648,
    H_ref = -934873,
    S_ref = 306.269,
    a1 = 6.391e-05,
    a2 = 11671.1,
    a3 = -7.2638e-05,
    a4 = -164519,
    c1 = 252.5931,
    c2 = 57409,
    w_ref = -42720,
    z = 0,
    a_0 = 3.72,
    MM = 0.13211) "GLUTARIC-ACID,aq";

  //197
  constant DataRecord GLYCINE_aq(
    R = Modelica.Constants.R/GLYCINE_aq.MM,
    G_ref = -370778,
    H_ref = -513988,
    S_ref = 158.323,
    a1 = 3.1818e-05,
    a2 = 2963.32,
    a3 = 0.00045655,
    a4 = -128520,
    c1 = 65.5838,
    c2 = -173402,
    w_ref = -97490,
    z = 0,
    a_0 = 3.72,
    MM = 0.075066) "GLYCINE,aq";

  //198
  constant DataRecord GLYCINATE_aq(
    R = Modelica.Constants.R/GLYCINATE_aq.MM,
    G_ref = -314972,
    H_ref = -469863,
    S_ref = 119.16,
    a1 = 3.4138e-05,
    a2 = 4924.4,
    a3 = 8.0375e-05,
    a4 = -136629,
    c1 = 54.1364,
    c2 = -175138,
    w_ref = 501030,
    z = -1,
    a_0 = 3.72,
    MM = 0.074058) "GLYCINATE,aq";

  //199
  constant DataRecord GLYCOLATE_aq(
    R = Modelica.Constants.R/GLYCOLATE_aq.MM,
    G_ref = -506975,
    H_ref = -647265,
    S_ref = 109.621,
    a1 = 3.1944e-05,
    a2 = 4626.25,
    a3 = 4.1145e-05,
    a4 = -135394,
    c1 = 108.9777,
    c2 = -168494,
    w_ref = 516050,
    z = -1,
    a_0 = 3.72,
    MM = 0.075044) "GLYCOLATE,aq";

  //200
  constant DataRecord GLYCOLIC_ACID_aq(
    R = Modelica.Constants.R/GLYCOLIC_ACID_aq.MM,
    G_ref = -528858,
    H_ref = -648060,
    S_ref = 180.33,
    a1 = 3.6587e-05,
    a2 = 5083.39,
    a3 = 0.00016829,
    a4 = -137285,
    c1 = 182.7023,
    c2 = -19477,
    w_ref = -126230,
    z = 0,
    a_0 = 3.72,
    MM = 0.076052) "GLYCOLIC-ACID,aq";

  //201
  constant DataRecord DIGLYCINE_aq(
    R = Modelica.Constants.R/DIGLYCINE_aq.MM,
    G_ref = -489612,
    H_ref = -734878,
    S_ref = 226.773,
    a1 = 5.0372e-05,
    a2 = 6671.26,
    a3 = 0.00048729,
    a4 = -143850,
    c1 = 152.0842,
    c2 = -57208,
    w_ref = -201170,
    z = 0,
    a_0 = 3.72,
    MM = 0.13212) "DIGLYCINE,aq";

  //202
  constant DataRecord Ga_p3(
    R = Modelica.Constants.R/Ga_p3.MM,
    G_ref = -158992,
    H_ref = -211710,
    S_ref = -330.536,
    a1 = -1.1873e-05,
    a2 = -6153.2,
    a3 = 0.00048198,
    a4 = -90830,
    c1 = 67.9298,
    c2 = -352820,
    w_ref = 1162610,
    z = 3,
    a_0 = 3.72,
    MM = 0.06972) "Ga+3";

  //203
  constant DataRecord GaOH_p2(
    R = Modelica.Constants.R/GaOH_p2.MM,
    G_ref = -381162,
    H_ref = -488273,
    S_ref = -281.583,
    a1 = 6.9337e-06,
    a2 = -1561.51,
    a3 = 0.0003017,
    a4 = -109813,
    c1 = 202.4512,
    c2 = 209681,
    w_ref = 866090,
    z = 2,
    a_0 = 3.72,
    MM = 0.086728) "GaOH+2";

  //204
  constant DataRecord GaO_p1(
    R = Modelica.Constants.R/GaO_p1.MM,
    G_ref = -362334,
    H_ref = -422166,
    S_ref = -122.173,
    a1 = 8.9115e-06,
    a2 = -1079.72,
    a3 = 0.00028296,
    a4 = -111805,
    c1 = -3.1393,
    c2 = -361343,
    w_ref = 417810,
    z = 1,
    a_0 = 3.72,
    MM = 0.08572) "GaO+";

  //205
  constant DataRecord HGaO2_aq(
    R = Modelica.Constants.R/HGaO2_aq.MM,
    G_ref = -574463,
    H_ref = -663582,
    S_ref = 12.134,
    a1 = 1.4422e-05,
    a2 = 264.76,
    a3 = 0.00023037,
    a4 = -117365,
    c1 = -50.926,
    c2 = -387761,
    w_ref = -18370,
    z = 0,
    a_0 = 3.72,
    MM = 0.10273) "HGaO2,aq";

  //206
  constant DataRecord GaO2_n1(
    R = Modelica.Constants.R/GaO2_n1.MM,
    G_ref = -538481,
    H_ref = -643918,
    S_ref = -41.84,
    a1 = 1.5166e-05,
    a2 = 448.57,
    a3 = 0.00022268,
    a4 = -118123,
    c1 = 60.2404,
    c2 = -245429,
    w_ref = 743790,
    z = -1,
    a_0 = 3.72,
    MM = 0.10172) "GaO2-";

  //207
  constant DataRecord Gd_p2(
    R = Modelica.Constants.R/Gd_p2.MM,
    G_ref = -295390,
    H_ref = -282002,
    S_ref = -17.991,
    a1 = 3.93e-08,
    a2 = -3243.27,
    a3 = 0.0003674,
    a4 = -102859,
    c1 = 41.8174,
    c2 = -221568,
    w_ref = 469280,
    z = 2,
    a_0 = 3.72,
    MM = 0.15725) "Gd+2";

  //208
  constant DataRecord Gd_p3(
    R = Modelica.Constants.R/Gd_p3.MM,
    G_ref = -663582,
    H_ref = -687013,
    S_ref = -205.853,
    a1 = -1.2456e-05,
    a2 = -6297.17,
    a3 = 0.00048809,
    a4 = -90241,
    c1 = 27.4496,
    c2 = -432935,
    w_ref = 973410,
    z = 3,
    a_0 = 3.72,
    MM = 0.15725) "Gd+3";

  //209
  constant DataRecord Gd_p4(
    R = Modelica.Constants.R/Gd_p4.MM,
    G_ref = 52718,
    H_ref = -23849,
    S_ref = -450.198,
    a1 = -1.8214e-05,
    a2 = -7700.74,
    a3 = 0.00054267,
    a4 = -84433,
    c1 = 166.688,
    c2 = -146566,
    w_ref = 1590460,
    z = 4,
    a_0 = 3.72,
    MM = 0.15725) "Gd+4";

  //210
  constant DataRecord H_p1(
    R = Modelica.Constants.R/H_p1.MM,
    G_ref = 0,
    H_ref = 0,
    S_ref = 0,
    a1 = 0,
    a2 = 0,
    a3 = 0,
    a4 = 0,
    c1 = 0,
    c2 = 0,
    w_ref = 0,
    z = 1,
    a_0 = 9,
    MM = 0.001008) "H+";

  //211
  constant DataRecord H_ADIPATE_aq(
    R = Modelica.Constants.R/H_ADIPATE_aq.MM,
    G_ref = -697766,
    H_ref = -950312,
    S_ref = 250.203,
    a1 = 6.6858e-05,
    a2 = 12340.12,
    a3 = -8.8061e-05,
    a4 = -167285,
    c1 = 210.1121,
    c2 = -153486,
    w_ref = 302670,
    z = -1,
    a_0 = 3.72,
    MM = 0.14513) "H-ADIPATE,aq";

  //212
  constant DataRecord H_AZELATE_aq(
    R = Modelica.Constants.R/H_AZELATE_aq.MM,
    G_ref = -658938,
    H_ref = -1008218,
    S_ref = 334.72,
    a1 = 9.3283e-05,
    a2 = 18329.14,
    a3 = -0.00022407,
    a4 = -192041,
    c1 = 458.2785,
    c2 = -121181,
    w_ref = 174640,
    z = -1,
    a_0 = 3.72,
    MM = 0.18721) "H-AZELATE,aq";

  //213
  constant DataRecord H_GLUTARATE_aq(
    R = Modelica.Constants.R/H_GLUTARATE_aq.MM,
    G_ref = -714878,
    H_ref = -935417,
    S_ref = 221.334,
    a1 = 5.7725e-05,
    a2 = 10270.46,
    a3 = -4.1099e-05,
    a4 = -158728,
    c1 = 148.8512,
    c2 = -161599,
    w_ref = 346390,
    z = -1,
    a_0 = 3.72,
    MM = 0.13111) "H-GLUTARATE,aq";

  //214
  constant DataRecord H_MALONATE_aq(
    R = Modelica.Constants.R/H_MALONATE_aq.MM,
    G_ref = -717682,
    H_ref = -869644,
    S_ref = 178.657,
    a1 = 4.1449e-05,
    a2 = 6581.14,
    a3 = 4.279e-05,
    a4 = -143474,
    c1 = 49.9331,
    c2 = -174632,
    w_ref = 411120,
    z = -1,
    a_0 = 3.72,
    MM = 0.10305) "H-MALONATE,aq";

  //215
  constant DataRecord H_OXALATE_aq(
    R = Modelica.Constants.R/H_OXALATE_aq.MM,
    G_ref = -698339,
    H_ref = -818390,
    S_ref = 150.624,
    a1 = 3.3169e-05,
    a2 = 4703.07,
    a3 = 8.5772e-05,
    a4 = -135712,
    c1 = 27.5788,
    c2 = -265885,
    w_ref = 453170,
    z = -1,
    a_0 = 3.72,
    MM = 0.089028) "H-OXALATE,aq";

  //216
  constant DataRecord H_PIMELATE_aq(
    R = Modelica.Constants.R/H_PIMELATE_aq.MM,
    G_ref = -695632,
    H_ref = -979223,
    S_ref = 282.42,
    a1 = 7.5631e-05,
    a2 = 14326.69,
    a3 = -0.00013281,
    a4 = -175494,
    c1 = 278.6728,
    c2 = -144407,
    w_ref = 253930,
    z = -1,
    a_0 = 3.72,
    MM = 0.15916) "H-PIMELATE,aq";

  //217
  constant DataRecord H_SEBACATE_aq(
    R = Modelica.Constants.R/H_SEBACATE_aq.MM,
    G_ref = -648687,
    H_ref = -1030226,
    S_ref = 362.753,
    a1 = 0.00010219,
    a2 = 20346.88,
    a3 = -0.00026992,
    a4 = -200384,
    c1 = 540.8795,
    c2 = -110432,
    w_ref = 132130,
    z = -1,
    a_0 = 3.72,
    MM = 0.20123) "H-SEBACATE,aq";

  //218
  constant DataRecord H_SUBERATE_aq(
    R = Modelica.Constants.R/H_SUBERATE_aq.MM,
    G_ref = -679231,
    H_ref = -996336,
    S_ref = 306.269,
    a1 = 8.4383e-05,
    a2 = 16310.74,
    a3 = -0.00017795,
    a4 = -183699,
    c1 = 375.737,
    c2 = -131934,
    w_ref = 217780,
    z = -1,
    a_0 = 3.72,
    MM = 0.17318) "H-SUBERATE,aq";

  //219
  constant DataRecord H_SUCCINATE_aq(
    R = Modelica.Constants.R/H_SUCCINATE_aq.MM,
    G_ref = -719899,
    H_ref = -909392,
    S_ref = 189.117,
    a1 = 4.8793e-05,
    a2 = 8244.86,
    a3 = 5.13e-06,
    a4 = -150352,
    c1 = 109.2735,
    c2 = -167075,
    w_ref = 395220,
    z = -1,
    a_0 = 3.72,
    MM = 0.11708) "H-SUCCINATE,aq";

  //220
  constant DataRecord H2_aq(
    R = Modelica.Constants.R/H2_aq.MM,
    G_ref = 17723,
    H_ref = -4184,
    S_ref = 57.739,
    a1 = 2.1517e-05,
    a2 = 1998.19,
    a3 = 0.00016204,
    a4 = -124533,
    c1 = 115.5834,
    c2 = 213091,
    w_ref = -87450,
    z = 0,
    a_0 = 3.72,
    MM = 0.0020158) "H2,aq";

  //221
  constant DataRecord H2AsO3_n1(
    R = Modelica.Constants.R/H2AsO3_n1.MM,
    G_ref = -587141,
    H_ref = -714795,
    S_ref = 110.458,
    a1 = 2.424e-05,
    a2 = 2662.95,
    a3 = 0.00013592,
    a4 = -127281,
    c1 = 66.1206,
    c2 = -151683,
    w_ref = 514840,
    z = -1,
    a_0 = 3.72,
    MM = 0.12494) "H2AsO3-";

  //222
  constant DataRecord H2AsO4_n1(
    R = Modelica.Constants.R/H2AsO4_n1.MM,
    G_ref = -753162,
    H_ref = -909560,
    S_ref = 117.152,
    a1 = 2.8212e-05,
    a2 = 3633.18,
    a3 = 9.7701e-05,
    a4 = -131290,
    c1 = 70.7958,
    c2 = -132076,
    w_ref = 504380,
    z = 1,
    a_0 = 4.25,
    MM = 0.14094) "H2AsO4-";

  //223
  constant DataRecord NaH2AsO4_aq(
    R = Modelica.Constants.R/NaH2AsO4_aq.MM,
    G_ref = -1004909,
    H_ref = -1140592,
    S_ref = 172.381,
    a1 = 3.0866e-05,
    a2 = 4282.32,
    a3 = 4.7785e-05,
    a4 = -133980,
    c1 = 122.172,
    c2 = 194589,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.16393) "NaH2AsO4,aq";

  //224
  constant DataRecord KH2AsO4_aq(
    R = Modelica.Constants.R/KH2AsO4_aq.MM,
    G_ref = -1024808,
    H_ref = -1147889,
    S_ref = 228.028,
    a1 = 3.7353e-05,
    a2 = 5866.8,
    a3 = -2.0251e-05,
    a4 = -140528,
    c1 = 81.6014,
    c2 = 58823,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.18004) "KH2AsO4,aq";

  //225
  constant DataRecord MgH2AsO4_p1(
    R = Modelica.Constants.R/MgH2AsO4_p1.MM,
    G_ref = -1217167,
    H_ref = -1392088,
    S_ref = -43.514,
    a1 = 1.8872e-05,
    a2 = 1352.9,
    a3 = 0.00017357,
    a4 = -121867,
    c1 = 126.5434,
    c2 = 304532,
    w_ref = 293260,
    z = 1,
    a_0 = 3.72,
    MM = 0.16525) "MgH2AsO4+";

  //226
  constant DataRecord CaH2AsO4_p1(
    R = Modelica.Constants.R/CaH2AsO4_p1.MM,
    G_ref = -1314487,
    H_ref = -1456061,
    S_ref = 77.404,
    a1 = 2.048e-05,
    a2 = 1745.73,
    a3 = 0.0001567,
    a4 = -123491,
    c1 = 130.8395,
    c2 = 262538,
    w_ref = 110420,
    z = 1,
    a_0 = 3.72,
    MM = 0.18102) "CaH2AsO4+";

  //227
  constant DataRecord SrH2AsO4_p1(
    R = Modelica.Constants.R/SrH2AsO4_p1.MM,
    G_ref = -1321709,
    H_ref = -1455367,
    S_ref = 114.223,
    a1 = 2.0706e-05,
    a2 = 1800.84,
    a3 = 0.00015434,
    a4 = -123721,
    c1 = 121.5565,
    c2 = 214212,
    w_ref = 54430,
    z = 1,
    a_0 = 3.72,
    MM = 0.22856) "SrH2AsO4+";

  //228
  constant DataRecord MnH2AsO4_p1(
    R = Modelica.Constants.R/MnH2AsO4_p1.MM,
    G_ref = -989441,
    H_ref = -1133144,
    S_ref = 60.668,
    a1 = 2.1178e-05,
    a2 = 1916.23,
    a3 = 0.00014938,
    a4 = -124198,
    c1 = 148.1621,
    c2 = 328314,
    w_ref = 135730,
    z = 1,
    a_0 = 3.72,
    MM = 0.19588) "MnH2AsO4+";

  //229
  constant DataRecord FeH2AsO4_p1(
    R = Modelica.Constants.R/FeH2AsO4_p1.MM,
    G_ref = -860620,
    H_ref = -1019720,
    S_ref = 4.184,
    a1 = 1.8213e-05,
    a2 = 1192.02,
    a3 = 0.00018047,
    a4 = -121202,
    c1 = 118.5256,
    c2 = 255442,
    w_ref = 221040,
    z = 1,
    a_0 = 3.72,
    MM = 0.19679) "FeH2AsO4+";

  //230
  constant DataRecord CoH2AsO4_p1(
    R = Modelica.Constants.R/CoH2AsO4_p1.MM,
    G_ref = -809135,
    H_ref = -972182,
    S_ref = -6.276,
    a1 = 1.6162e-05,
    a2 = 691.03,
    a3 = 0.00020198,
    a4 = -119131,
    c1 = 117.6298,
    c2 = 257362,
    w_ref = 236980,
    z = 1,
    a_0 = 3.72,
    MM = 0.19987) "CoH2AsO4+";

  //231
  constant DataRecord NiH2AsO4_p1(
    R = Modelica.Constants.R/NiH2AsO4_p1.MM,
    G_ref = -808123,
    H_ref = -978156,
    S_ref = -29.706,
    a1 = 1.405e-05,
    a2 = 175.1,
    a3 = 0.00022413,
    a4 = -116997,
    c1 = 91.999,
    c2 = 182573,
    w_ref = 272630,
    z = 1,
    a_0 = 3.72,
    MM = 0.19965) "NiH2AsO4+";

  //232
  constant DataRecord CuH2AsO4_p1(
    R = Modelica.Constants.R/CuH2AsO4_p1.MM,
    G_ref = -698167,
    H_ref = -855134,
    S_ref = 17.154,
    a1 = 1.6616e-05,
    a2 = 801.82,
    a3 = 0.00019723,
    a4 = -119591,
    c1 = 132.9462,
    c2 = 297629,
    w_ref = 201380,
    z = 1,
    a_0 = 3.72,
    MM = 0.20448) "CuH2AsO4+";

  //233
  constant DataRecord ZnH2AsO4_p1(
    R = Modelica.Constants.R/ZnH2AsO4_p1.MM,
    G_ref = -903439,
    H_ref = -1068476,
    S_ref = -1.255,
    a1 = 1.6902e-05,
    a2 = 871.78,
    a3 = 0.00019422,
    a4 = -119880,
    c1 = 132.4751,
    c2 = 304725,
    w_ref = 229490,
    z = 1,
    a_0 = 3.72,
    MM = 0.20631) "ZnH2AsO4+";

  //234
  constant DataRecord PbH2AsO4_p1(
    R = Modelica.Constants.R/PbH2AsO4_p1.MM,
    G_ref = -786157,
    H_ref = -901974,
    S_ref = 187.025,
    a1 = 2.1489e-05,
    a2 = 1992.04,
    a3 = 0.00014612,
    a4 = -124512,
    c1 = 116.5047,
    c2 = 163398,
    w_ref = -55560,
    z = 1,
    a_0 = 3.72,
    MM = 0.34813) "PbH2AsO4+";

  //235
  constant DataRecord AlH2AsO4_p2(
    R = Modelica.Constants.R/AlH2AsO4_p2.MM,
    G_ref = -1255108,
    H_ref = -1467400,
    S_ref = -238.488,
    a1 = 6.063e-06,
    a2 = -1775.69,
    a3 = 0.00030789,
    a4 = -108935,
    c1 = 38.1003,
    c2 = 169787,
    w_ref = 816210,
    z = 2,
    a_0 = 3.72,
    MM = 0.16792) "AlH2AsO4+2";

  //236
  constant DataRecord FeH2AsO4_p2(
    R = Modelica.Constants.R/FeH2AsO4_p2.MM,
    G_ref = -794747,
    H_ref = -985478,
    S_ref = -166.942,
    a1 = 1.0419e-05,
    a2 = -711.7,
    a3 = 0.00026221,
    a4 = -113332,
    c1 = 127.6994,
    c2 = 436308,
    w_ref = 708230,
    z = 2,
    a_0 = 3.72,
    MM = 0.19679) "FeH2AsO4+2";

  //237
  constant DataRecord NaHAsO4_n1(
    R = Modelica.Constants.R/NaHAsO4_n1.MM,
    G_ref = -979165,
    H_ref = -1141973,
    S_ref = 81.17,
    a1 = 1.8677e-05,
    a2 = 1305.41,
    a3 = 0.0001756,
    a4 = -121671,
    c1 = 1.3744,
    c2 = -33531,
    w_ref = 555380,
    z = 1,
    a_0 = 3.72,
    MM = 0.16292) "NaHAsO4-";

  //238
  constant DataRecord KHAsO4_n1(
    R = Modelica.Constants.R/KHAsO4_n1.MM,
    G_ref = -998947,
    H_ref = -1151633,
    S_ref = 128.449,
    a1 = 2.3967e-05,
    a2 = 2597.34,
    a3 = 0.00012014,
    a4 = -127014,
    c1 = -10.055,
    c2 = -93872,
    w_ref = 483710,
    z = 1,
    a_0 = 3.72,
    MM = 0.17903) "KHAsO4-";

  //239
  constant DataRecord MgHAsO4_aq(
    R = Modelica.Constants.R/MgHAsO4_aq.MM,
    G_ref = -1182629,
    H_ref = -1365892,
    S_ref = -71.546,
    a1 = 7.0647e-06,
    a2 = -1530.93,
    a3 = 0.00029738,
    a4 = -109947,
    c1 = 17.3305,
    c2 = -156231,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.16424) "MgHAsO4,aq";

  //240
  constant DataRecord CaHAsO4_aq(
    R = Modelica.Constants.R/CaHAsO4_aq.MM,
    G_ref = -1280463,
    H_ref = -1440723,
    S_ref = 14.644,
    a1 = 8.4634e-06,
    a2 = -1189.09,
    a3 = 0.00028271,
    a4 = -111357,
    c1 = 11.7529,
    c2 = -174891,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.18001) "CaHAsO4,aq";

  //241
  constant DataRecord SrHAsO4_aq(
    R = Modelica.Constants.R/SrHAsO4_aq.MM,
    G_ref = -1287626,
    H_ref = -1443133,
    S_ref = 41.003,
    a1 = 8.7241e-06,
    a2 = -1125.5,
    a3 = 0.00027998,
    a4 = -111621,
    c1 = 5.335,
    c2 = -196397,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.22755) "SrHAsO4,aq";

  //242
  constant DataRecord MnHAsO4_aq(
    R = Modelica.Constants.R/MnHAsO4_aq.MM,
    G_ref = -960496,
    H_ref = -1121450,
    S_ref = 2.929,
    a1 = 8.8483e-06,
    a2 = -1095.37,
    a3 = 0.00027868,
    a4 = -111746,
    c1 = 20.4886,
    c2 = -145687,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.19487) "MnHAsO4,aq";

  //243
  constant DataRecord FeHAsO4_aq(
    R = Modelica.Constants.R/FeHAsO4_aq.MM,
    G_ref = -824085,
    H_ref = -995616,
    S_ref = -37.238,
    a1 = 6.8044e-06,
    a2 = -1594.52,
    a3 = 0.00030011,
    a4 = -109684,
    c1 = 10.8106,
    c2 = -178071,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.19578) "FeHAsO4,aq";

  //244
  constant DataRecord CoHAsO4_aq(
    R = Modelica.Constants.R/CoHAsO4_aq.MM,
    G_ref = -784567,
    H_ref = -959140,
    S_ref = -44.769,
    a1 = 5.4819e-06,
    a2 = -1917.53,
    a3 = 0.00031398,
    a4 = -108349,
    c1 = 11.0654,
    c2 = -177192,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.19886) "CoHAsO4,aq";

  //245
  constant DataRecord NiHAsO4_aq(
    R = Modelica.Constants.R/NiHAsO4_aq.MM,
    G_ref = -774387,
    H_ref = -953931,
    S_ref = -61.505,
    a1 = 4.079e-06,
    a2 = -2260.2,
    a3 = 0.00032869,
    a4 = -106930,
    c1 = 1.1326,
    c2 = -210455,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.19864) "NiHAsO4,aq";

  //246
  constant DataRecord CuHAsO4_aq(
    R = Modelica.Constants.R/CuHAsO4_aq.MM,
    G_ref = -669624,
    H_ref = -840135,
    S_ref = -28.033,
    a1 = 5.8425e-06,
    a2 = -1829.24,
    a3 = 0.0003102,
    a4 = -108713,
    c1 = 16.4134,
    c2 = -159327,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.20347) "CuHAsO4,aq";

  //247
  constant DataRecord ZnHAsO4_aq(
    R = Modelica.Constants.R/ZnHAsO4_aq.MM,
    G_ref = -877916,
    H_ref = -1054904,
    S_ref = -41.422,
    a1 = 5.9626e-06,
    a2 = -1799.96,
    a3 = 0.00030894,
    a4 = -108834,
    c1 = 17.3561,
    c2 = -156147,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.2053) "ZnHAsO4,aq";

  //248
  constant DataRecord PbHAsO4_aq(
    R = Modelica.Constants.R/PbHAsO4_aq.MM,
    G_ref = -753618,
    H_ref = -897506,
    S_ref = 92.885,
    a1 = 9.4496e-06,
    a2 = -948.51,
    a3 = 0.00027237,
    a4 = -112353,
    c1 = -1.4142,
    c2 = -218949,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.34712) "PbHAsO4,aq";

  //249
  constant DataRecord AlHAsO4_p1(
    R = Modelica.Constants.R/AlHAsO4_p1.MM,
    G_ref = -1237979,
    H_ref = -1436694,
    S_ref = -192.882,
    a1 = 4.4233e-06,
    a2 = -2176.1,
    a3 = 0.00032508,
    a4 = -107278,
    c1 = -118.4072,
    c2 = -445387,
    w_ref = 519650,
    z = 1,
    a_0 = 3.72,
    MM = 0.16691) "AlHAsO4+";

  //250
  constant DataRecord FeHAsO4_p1(
    R = Modelica.Constants.R/FeHAsO4_p1.MM,
    G_ref = -787379,
    H_ref = -971571,
    S_ref = -145.185,
    a1 = 6.0856e-06,
    a2 = -1769.83,
    a3 = 0.00030765,
    a4 = -108960,
    c1 = -76.3496,
    c2 = -326896,
    w_ref = 447440,
    z = 1,
    a_0 = 3.72,
    MM = 0.19578) "FeHAsO4+";

  //251
  constant DataRecord NaAsO4_n2(
    R = Modelica.Constants.R/NaAsO4_n2.MM,
    G_ref = -935961,
    H_ref = -1062498,
    S_ref = 202.924,
    a1 = 1.1816e-06,
    a2 = -2967.71,
    a3 = 0.00035908,
    a4 = -104006,
    c1 = -170.1131,
    c2 = -455052,
    w_ref = 1049560,
    z = 2,
    a_0 = 3.72,
    MM = 0.16191) "NaAsO4-2";

  //252
  constant DataRecord KAsO4_n2(
    R = Modelica.Constants.R/KAsO4_n2.MM,
    G_ref = -955743,
    H_ref = -1042021,
    S_ref = 351.456,
    a1 = 4.653e-06,
    a2 = -2119.61,
    a3 = 0.00032267,
    a4 = -107512,
    c1 = -167.4437,
    c2 = -515385,
    w_ref = 824830,
    z = 2,
    a_0 = 3.72,
    MM = 0.17802) "KAsO4-2";

  //253
  constant DataRecord MgAsO4_n1(
    R = Modelica.Constants.R/MgAsO4_n1.MM,
    G_ref = -1135847,
    H_ref = -1276819,
    S_ref = 70.71,
    a1 = -5.558e-06,
    a2 = -4613.7,
    a3 = 0.00042975,
    a4 = -97203,
    c1 = -168.1717,
    c2 = -595885,
    w_ref = 571620,
    z = 1,
    a_0 = 3.72,
    MM = 0.16323) "MgAsO4-";

  //254
  constant DataRecord CaAsO4_n1(
    R = Modelica.Constants.R/CaAsO4_n1.MM,
    G_ref = -1233916,
    H_ref = -1351750,
    S_ref = 156.9,
    a1 = -5.7237e-06,
    a2 = -4652.61,
    a3 = 0.0004315,
    a4 = -97027,
    c1 = -161.6781,
    c2 = -614546,
    w_ref = 440620,
    z = 1,
    a_0 = 3.72,
    MM = 0.179) "CaAsO4-";

  //255
  constant DataRecord SrAsO4_n1(
    R = Modelica.Constants.R/SrAsO4_n1.MM,
    G_ref = -1239481,
    H_ref = -1352432,
    S_ref = 183.678,
    a1 = -5.8074e-06,
    a2 = -4673.53,
    a3 = 0.00043238,
    a4 = -96943,
    c1 = -164.3601,
    c2 = -636010,
    w_ref = 400070,
    z = 1,
    a_0 = 3.72,
    MM = 0.22654) "SrAsO4-";

  //256
  constant DataRecord MnAsO4_n1(
    R = Modelica.Constants.R/MnAsO4_n1.MM,
    G_ref = -913321,
    H_ref = -1051845,
    S_ref = 78.241,
    a1 = -5.2472e-06,
    a2 = -4537.97,
    a3 = 0.00042649,
    a4 = -97516,
    c1 = -163.9626,
    c2 = -585300,
    w_ref = 560240,
    z = 1,
    a_0 = 3.72,
    MM = 0.19386) "MnAsO4-";

  //257
  constant DataRecord FeAsO4_n1(
    R = Modelica.Constants.R/FeAsO4_n1.MM,
    G_ref = -781019,
    H_ref = -914376,
    S_ref = 90.793,
    a1 = -5.7124e-06,
    a2 = -4651.77,
    a3 = 0.00043138,
    a4 = -97044,
    c1 = -171.8787,
    c2 = -616554,
    w_ref = 541070,
    z = 1,
    a_0 = 3.72,
    MM = 0.19477) "FeAsO4-";

  //258
  constant DataRecord CoAsO4_n1(
    R = Modelica.Constants.R/CoAsO4_n1.MM,
    G_ref = -741363,
    H_ref = -880653,
    S_ref = 73.638,
    a1 = -5.8881e-06,
    a2 = -4694.45,
    a3 = 0.00043322,
    a4 = -96868,
    c1 = -174.0209,
    c2 = -616847,
    w_ref = 567140,
    z = 1,
    a_0 = 3.72,
    MM = 0.19785) "CoAsO4-";

  //259
  constant DataRecord NiAsO4_n1(
    R = Modelica.Constants.R/NiAsO4_n1.MM,
    G_ref = -737656,
    H_ref = -884163,
    S_ref = 49.371,
    a1 = -6.0383e-06,
    a2 = -4731.27,
    a3 = 0.00043479,
    a4 = -96717,
    c1 = -187.347,
    c2 = -650110,
    w_ref = 603920,
    z = 1,
    a_0 = 3.72,
    MM = 0.19763) "NiAsO4-";

  //260
  constant DataRecord CuAsO4_n1(
    R = Modelica.Constants.R/CuAsO4_n1.MM,
    G_ref = -634893,
    H_ref = -768212,
    S_ref = 96.65,
    a1 = -5.9321e-06,
    a2 = -4705.33,
    a3 = 0.00043368,
    a4 = -96822,
    c1 = -165.4395,
    c2 = -598940,
    w_ref = 532040,
    z = 1,
    a_0 = 3.72,
    MM = 0.20246) "CuAsO4-";

  //261
  constant DataRecord ZnAsO4_n1(
    R = Modelica.Constants.R/ZnAsO4_n1.MM,
    G_ref = -837306,
    H_ref = -978441,
    S_ref = 79.078,
    a1 = -5.8187e-06,
    a2 = -4677.29,
    a3 = 0.00043248,
    a4 = -96939,
    c1 = -166.9709,
    c2 = -595802,
    w_ref = 558900,
    z = 1,
    a_0 = 3.72,
    MM = 0.20429) "ZnAsO4-";

  //262
  constant DataRecord PbAsO4_n1(
    R = Modelica.Constants.R/PbAsO4_n1.MM,
    G_ref = -710414,
    H_ref = -813524,
    S_ref = 229.702,
    a1 = -5.9028e-06,
    a2 = -4698.21,
    a3 = 0.00043337,
    a4 = -96851,
    c1 = -164.7199,
    c2 = -658603,
    w_ref = 330700,
    z = 1,
    a_0 = 3.72,
    MM = 0.34611) "PbAsO4-";

  //263
  constant DataRecord AlAsO4_aq(
    R = Modelica.Constants.R/AlAsO4_aq.MM,
    G_ref = -1194775,
    H_ref = -1382385,
    S_ref = -155.645,
    a1 = 1.7489e-06,
    a2 = -2829.22,
    a3 = 0.00035313,
    a4 = -104579,
    c1 = -127.4237,
    c2 = -640654,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.1659) "AlAsO4,aq";

  //264
  constant DataRecord FeAsO4_aq(
    R = Modelica.Constants.R/FeAsO4_aq.MM,
    G_ref = -744175,
    H_ref = -915894,
    S_ref = -103.345,
    a1 = -1.636e-07,
    a2 = -3296.16,
    a3 = 0.00037318,
    a4 = -102646,
    c1 = -190.0749,
    c2 = -850314,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.19477) "FeAsO4,aq";

  //265
  constant DataRecord NaH2AsO3_aq(
    R = Modelica.Constants.R/NaH2AsO3_aq.MM,
    G_ref = -850482,
    H_ref = -957433,
    S_ref = 166.105,
    a1 = 2.604e-05,
    a2 = 3103.27,
    a3 = 9.8399e-05,
    a4 = -129106,
    c1 = 98.2571,
    c2 = 114558,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.14793) "NaH2AsO3,aq";

  //266
  constant DataRecord AgH2AsO3_aq(
    R = Modelica.Constants.R/AgH2AsO3_aq.MM,
    G_ref = -516833,
    H_ref = -615345,
    S_ref = 185.351,
    a1 = 2.6237e-05,
    a2 = 3151.81,
    a3 = 9.6324e-05,
    a4 = -129307,
    c1 = 91.6087,
    c2 = 92299,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.23281) "AgH2AsO3,aq";

  //267
  constant DataRecord MgH2AsO3_p1(
    R = Modelica.Constants.R/MgH2AsO3_p1.MM,
    G_ref = -1051899,
    H_ref = -1197390,
    S_ref = -47.279,
    a1 = 1.4065e-05,
    a2 = 178.66,
    a3 = 0.00022397,
    a4 = -117014,
    c1 = 102.0854,
    c2 = 224472,
    w_ref = 299110,
    z = 1,
    a_0 = 3.72,
    MM = 0.14925) "MgH2AsO3+";

  //268
  constant DataRecord CaH2AsO3_p1(
    R = Modelica.Constants.R/CaH2AsO3_p1.MM,
    G_ref = -1150274,
    H_ref = -1263032,
    S_ref = 71.546,
    a1 = 1.5684e-05,
    a2 = 574.04,
    a3 = 0.000207,
    a4 = -118650,
    c1 = 106.0979,
    c2 = 182506,
    w_ref = 119370,
    z = 1,
    a_0 = 3.72,
    MM = 0.16502) "CaH2AsO3+";

  //269
  constant DataRecord SrH2AsO3_p1(
    R = Modelica.Constants.R/SrH2AsO3_p1.MM,
    G_ref = -1153115,
    H_ref = -1258137,
    S_ref = 107.947,
    a1 = 1.5913e-05,
    a2 = 630.11,
    a3 = 0.0002046,
    a4 = -118880,
    c1 = 96.7257,
    c2 = 134181,
    w_ref = 64350,
    z = 1,
    a_0 = 3.72,
    MM = 0.21256) "SrH2AsO3+";

  //270
  constant DataRecord BaH2AsO3_p1(
    R = Modelica.Constants.R/BaH2AsO3_p1.MM,
    G_ref = -1156102,
    H_ref = -1246564,
    S_ref = 167.778,
    a1 = 1.8675e-05,
    a2 = 1304.57,
    a3 = 0.00017563,
    a4 = -121671,
    c1 = 92.1819,
    c2 = 91002,
    w_ref = -26280,
    z = 1,
    a_0 = 3.72,
    MM = 0.26228) "BaH2AsO3+";

  //271
  constant DataRecord CuH2AsO3_p1(
    R = Modelica.Constants.R/CuH2AsO3_p1.MM,
    G_ref = -562196,
    H_ref = -690042,
    S_ref = 12.552,
    a1 = 1.1814e-05,
    a2 = -370.7,
    a3 = 0.00024758,
    a4 = -114742,
    c1 = 108.3447,
    c2 = 217568,
    w_ref = 208780,
    z = 1,
    a_0 = 3.72,
    MM = 0.18848) "CuH2AsO3+";

  //272
  constant DataRecord PbH2AsO3_p1(
    R = Modelica.Constants.R/PbH2AsO3_p1.MM,
    G_ref = -640729,
    H_ref = -728280,
    S_ref = 179.075,
    a1 = 1.6702e-05,
    a2 = 822.57,
    a3 = 0.00019632,
    a4 = -119675,
    c1 = 91.5041,
    c2 = 83345,
    w_ref = -43810,
    z = 1,
    a_0 = 3.72,
    MM = 0.33213) "PbH2AsO3+";

  //273
  constant DataRecord AlH2AsO3_p2(
    R = Modelica.Constants.R/AlH2AsO3_p2.MM,
    G_ref = -1115500,
    H_ref = -1295417,
    S_ref = -232.212,
    a1 = 1.2058e-06,
    a2 = -2961.85,
    a3 = 0.00035882,
    a4 = -104031,
    c1 = 15.0243,
    c2 = 89747,
    w_ref = 807090,
    z = 2,
    a_0 = 3.72,
    MM = 0.15192) "AlH2AsO3+2";

  //274
  constant DataRecord FeH2AsO3_p2(
    R = Modelica.Constants.R/FeH2AsO3_p2.MM,
    G_ref = -645976,
    H_ref = -805056,
    S_ref = -163.594,
    a1 = 5.5735e-06,
    a2 = -1894.93,
    a3 = 0.00031302,
    a4 = -108441,
    c1 = 104.282,
    c2 = 356309,
    w_ref = 702790,
    z = 2,
    a_0 = 3.72,
    MM = 0.18079) "FeH2AsO3+2";

  //275
  constant DataRecord H2P2O7_n2(
    R = Modelica.Constants.R/H2P2O7_n2.MM,
    G_ref = -2009994,
    H_ref = -2278606,
    S_ref = 163.176,
    a1 = 3.8059e-05,
    a2 = 6037.47,
    a3 = 3.18e-06,
    a4 = -141227,
    c1 = 170.0206,
    c2 = 23033,
    w_ref = 1096960,
    z = -2,
    a_0 = 3.72,
    MM = 0.17596) "H2P2O7-2";

  //276
  constant DataRecord H2PO4_n1(
    R = Modelica.Constants.R/H2PO4_n1.MM,
    G_ref = -1130266,
    H_ref = -1296287,
    S_ref = 90.374,
    a1 = 2.7144e-05,
    a2 = 3372.05,
    a3 = 0.00010804,
    a4 = -130214,
    c1 = 58.758,
    c2 = -186627,
    w_ref = 544050,
    z = -1,
    a_0 = 4.25,
    MM = 0.096986) "H2PO4-";

  //277
  constant DataRecord H2S_aq(
    R = Modelica.Constants.R/H2S_aq.MM,
    G_ref = -27920,
    H_ref = -37660,
    S_ref = 125.52,
    a1 = 2.7237e-05,
    a2 = 2833.57,
    a3 = 0.00024956,
    a4 = -127989,
    c1 = 135.1432,
    c2 = 197903,
    w_ref = -41840,
    z = 0,
    a_0 = 3.72,
    MM = 0.034076) "H2S,aq";

  //278
  constant DataRecord H3VO4_aq(
    R = Modelica.Constants.R/H3VO4_aq.MM,
    G_ref = -1042653,
    H_ref = -1190766,
    S_ref = 155.226,
    a1 = 2.8999e-05,
    a2 = 3823.84,
    a3 = 9.0542e-05,
    a4 = -132076,
    c1 = 64.0847,
    c2 = 35819,
    w_ref = -92800,
    z = 0,
    a_0 = 3.72,
    MM = 0.11796) "H3VO4,aq";

  //279
  constant DataRecord H2VO4_n1(
    R = Modelica.Constants.R/H2VO4_n1.MM,
    G_ref = -1020854,
    H_ref = -1174114,
    S_ref = 121.336,
    a1 = 2.6988e-05,
    a2 = 3332.68,
    a3 = 0.00010984,
    a4 = -130047,
    c1 = 159.4472,
    c2 = 178151,
    w_ref = 497810,
    z = -1,
    a_0 = 3.72,
    MM = 0.11696) "H2VO4-";

  //280
  constant DataRecord H3P2O7_n1(
    R = Modelica.Constants.R/H3P2O7_n1.MM,
    G_ref = -2023382,
    H_ref = -2276514,
    S_ref = 213.384,
    a1 = 3.8197e-05,
    a2 = 6071.9,
    a3 = 1.665e-06,
    a4 = -141369,
    c1 = 150.534,
    c2 = 191786,
    w_ref = 358490,
    z = -1,
    a_0 = 3.72,
    MM = 0.17696) "H3P2O7-";

  //281
  constant DataRecord H3PO4_aq(
    R = Modelica.Constants.R/H3PO4_aq.MM,
    G_ref = -1142650,
    H_ref = -1288337,
    S_ref = 158.992,
    a1 = 3.4613e-05,
    a2 = 5195.77,
    a3 = 3.6363e-05,
    a4 = -137754,
    c1 = 75.1898,
    c2 = 74170,
    w_ref = -92050,
    z = 0,
    a_0 = 3.72,
    MM = 0.097994) "H3PO4,aq";

  //282
  constant DataRecord HAsO4_n2(
    R = Modelica.Constants.R/HAsO4_n2.MM,
    G_ref = -714585,
    H_ref = -906338,
    S_ref = -1.674,
    a1 = 1.8407e-05,
    a2 = 1238.92,
    a3 = 0.00019185,
    a4 = -121390,
    c1 = 33.4335,
    c2 = -531795,
    w_ref = 1347120,
    z = -2,
    a_0 = 3.72,
    MM = 0.13993) "HAsO4-2";

  //283
  constant DataRecord HCO3_n1(
    R = Modelica.Constants.R/HCO3_n1.MM,
    G_ref = -586940,
    H_ref = -689933,
    S_ref = 98.45,
    a1 = 3.164e-05,
    a2 = 481.37,
    a3 = 5.1656e-05,
    a4 = -118265,
    c1 = 54.1389,
    c2 = -199071,
    w_ref = 532750,
    z = -1,
    a_0 = 4.25,
    MM = 0.061018) "HCO3-";

  //284
  constant DataRecord HCrO4_n1(
    R = Modelica.Constants.R/HCrO4_n1.MM,
    G_ref = -768601,
    H_ref = -878640,
    S_ref = 194.974,
    a1 = 3.4397e-05,
    a2 = 5143.18,
    a3 = 3.8384e-05,
    a4 = -137532,
    c1 = 65.0579,
    c2 = -114181,
    w_ref = 386180,
    z = -1,
    a_0 = 3.72,
    MM = 0.11701) "HCrO4-";

  //285
  constant DataRecord HEPTANE_aq(
    R = Modelica.Constants.R/HEPTANE_aq.MM,
    G_ref = 27070,
    H_ref = -221543,
    S_ref = 251.04,
    a1 = 8.0644e-05,
    a2 = 16435.17,
    a3 = -0.00040534,
    a4 = -184213,
    c1 = 472.3652,
    c2 = 1546913,
    w_ref = -380160,
    z = 0,
    a_0 = 3.72,
    MM = 0.1002) "HEPTANE,aq";

  //286
  constant DataRecord HEPTANOATE_aq(
    R = Modelica.Constants.R/HEPTANOATE_aq.MM,
    G_ref = -328549,
    H_ref = -609274,
    S_ref = 217.568,
    a1 = 7.551e-05,
    a2 = 13879.42,
    a3 = -2.5556e-05,
    a4 = -173649,
    c1 = 521.1553,
    c2 = -115403,
    w_ref = 352170,
    z = -1,
    a_0 = 3.72,
    MM = 0.12917) "HEPTANOATE,aq";

  //287
  constant DataRecord HEPTANOIC_ACID_aq(
    R = Modelica.Constants.R/HEPTANOIC_ACID_aq.MM,
    G_ref = -356477,
    H_ref = -607015,
    S_ref = 318.821,
    a1 = 8.2615e-05,
    a2 = 15909.7,
    a3 = -0.00016879,
    a4 = -182042,
    c1 = 527.8999,
    c2 = 396794,
    w_ref = -34390,
    z = 0,
    a_0 = 3.72,
    MM = 0.13018) "HEPTANOIC-ACID,aq";

  //288
  constant DataRecord HEXANE_aq(
    R = Modelica.Constants.R/HEXANE_aq.MM,
    G_ref = 18493,
    H_ref = -198322,
    S_ref = 221.334,
    a1 = 7.175e-05,
    a2 = 14264.39,
    a3 = -0.0003202,
    a4 = -175238,
    c1 = 424.5262,
    c2 = 1366227,
    w_ref = -335180,
    z = 0,
    a_0 = 3.72,
    MM = 0.086171) "HEXANE,aq";

  //289
  constant DataRecord HEXANOATE_aq(
    R = Modelica.Constants.R/HEXANOATE_aq.MM,
    G_ref = -336770,
    H_ref = -585216,
    S_ref = 189.535,
    a1 = 6.7237e-05,
    a2 = 12426.27,
    a3 = -9.0082e-05,
    a4 = -167640,
    c1 = 438.5313,
    c2 = -126152,
    w_ref = 394430,
    z = -1,
    a_0 = 3.72,
    MM = 0.11515) "HEXANOATE,aq";

  //290
  constant DataRecord HEXANOIC_ACID_aq(
    R = Modelica.Constants.R/HEXANOIC_ACID_aq.MM,
    G_ref = -364510,
    H_ref = -582789,
    S_ref = 290.788,
    a1 = 7.3935e-05,
    a2 = 13943.22,
    a3 = -0.00012427,
    a4 = -173912,
    c1 = 455.2958,
    c2 = 309156,
    w_ref = -52970,
    z = 0,
    a_0 = 3.72,
    MM = 0.11615) "HEXANOIC-ACID,aq";

  //291
  constant DataRecord HF_aq(
    R = Modelica.Constants.R/HF_aq.MM,
    G_ref = -299834,
    H_ref = -321478,
    S_ref = 94.14,
    a1 = 1.4541e-05,
    a2 = 294.64,
    a3 = 0.000229,
    a4 = -117491,
    c1 = 60.1019,
    c2 = -7648,
    w_ref = -290,
    z = 0,
    a_0 = 3.72,
    MM = 0.020008) "HF,aq";

  //292
  constant DataRecord HF2_n1(
    R = Modelica.Constants.R/HF2_n1.MM,
    G_ref = -578061,
    H_ref = -649943,
    S_ref = 92.466,
    a1 = 2.1867e-05,
    a2 = 2083.51,
    a3 = 0.00015869,
    a4 = -124888,
    c1 = -5.7534,
    c2 = -409923,
    w_ref = 541160,
    z = -1,
    a_0 = 3.72,
    MM = 0.039008) "HF2-";

  //293
  constant DataRecord HNO3_aq(
    R = Modelica.Constants.R/HNO3_aq.MM,
    G_ref = -103470,
    H_ref = -189995,
    S_ref = 178.657,
    a1 = 2.9967e-05,
    a2 = 4061.12,
    a3 = 8.1032e-05,
    a4 = -133060,
    c1 = 58.1187,
    c2 = 26443,
    w_ref = -128280,
    z = 0,
    a_0 = 3.72,
    MM = 0.063015) "HNO3,aq";

  //294
  constant DataRecord HO2_n1(
    R = Modelica.Constants.R/HO2_n1.MM,
    G_ref = -67321,
    H_ref = -160331,
    S_ref = 23.849,
    a1 = 1.1115e-05,
    a2 = -542.29,
    a3 = 0.000262,
    a4 = -114027,
    c1 = -6.3555,
    c2 = -445717,
    w_ref = 646390,
    z = -1,
    a_0 = 3.72,
    MM = 0.033008) "HO2-";

  //295
  constant DataRecord HPO4_n2(
    R = Modelica.Constants.R/HPO4_n2.MM,
    G_ref = -1089137,
    H_ref = -1292082,
    S_ref = -33.472,
    a1 = 1.5194e-05,
    a2 = 454.26,
    a3 = 0.00022273,
    a4 = -118152,
    c1 = 11.4462,
    c2 = -623847,
    w_ref = 1395910,
    z = -2,
    a_0 = 4,
    MM = 0.095978) "HPO4-2";

  //296
  constant DataRecord HS_n1(
    R = Modelica.Constants.R/HS_n1.MM,
    G_ref = 11966,
    H_ref = -16108,
    S_ref = 68.199,
    a1 = 2.097e-05,
    a2 = 2083.59,
    a3 = 0.00014546,
    a4 = -124888,
    c1 = 14.3093,
    c2 = -262337,
    w_ref = 602910,
    z = -1,
    a_0 = 3.5,
    MM = 0.033068) "HS-";

  //297
  constant DataRecord HSO3_n1(
    R = Modelica.Constants.R/HSO3_n1.MM,
    G_ref = -527728,
    H_ref = -626219,
    S_ref = 139.746,
    a1 = 2.8039e-05,
    a2 = 3590.54,
    a3 = 9.9458e-05,
    a4 = -131118,
    c1 = 65.6675,
    c2 = -138900,
    w_ref = 469990,
    z = -1,
    a_0 = 4.25,
    MM = 0.081068) "HSO3-";

  //298
  constant DataRecord HSO4_n1(
    R = Modelica.Constants.R/HSO4_n1.MM,
    G_ref = -755756,
    H_ref = -889100,
    S_ref = 125.52,
    a1 = 2.9199e-05,
    a2 = 3873.97,
    a3 = 8.8316e-05,
    a4 = -132290,
    c1 = 84.0821,
    c2 = -81797,
    w_ref = 491540,
    z = -1,
    a_0 = 3.72,
    MM = 0.097068) "HSO4-";

  //299
  constant DataRecord HSO5_n1(
    R = Modelica.Constants.R/HSO5_n1.MM,
    G_ref = -637516,
    H_ref = -775630,
    S_ref = 212.129,
    a1 = 3.7401e-05,
    a2 = 5875.59,
    a3 = 9.828e-06,
    a4 = -140557,
    c1 = 149.4751,
    c2 = 187523,
    w_ref = 360280,
    z = -1,
    a_0 = 3.72,
    MM = 0.11307) "HSO5-";

  //300
  constant DataRecord HSe_n1(
    R = Modelica.Constants.R/HSe_n1.MM,
    G_ref = 43932,
    H_ref = 15899,
    S_ref = 79.496,
    a1 = 2.0846e-05,
    a2 = 1834.89,
    a3 = 0.00016831,
    a4 = -123855,
    c1 = 28.9311,
    c2 = -295717,
    w_ref = 560660,
    z = -1,
    a_0 = 3.72,
    MM = 0.079968) "HSe-";

  //301
  constant DataRecord H2SeO3_aq(
    R = Modelica.Constants.R/H2SeO3_aq.MM,
    G_ref = -426140,
    H_ref = -507477,
    S_ref = 207.945,
    a1 = 3.3252e-05,
    a2 = 4863.48,
    a3 = 4.9451e-05,
    a4 = -136373,
    c1 = 96.2086,
    c2 = 173034,
    w_ref = -172670,
    z = 0,
    a_0 = 3.72,
    MM = 0.12898) "H2SeO3,aq";

  //302
  constant DataRecord HSeO3_n1(
    R = Modelica.Constants.R/HSeO3_n1.MM,
    G_ref = -411455,
    H_ref = -514548,
    S_ref = 135.143,
    a1 = 2.6803e-05,
    a2 = 3289.63,
    a3 = 0.00011109,
    a4 = -129867,
    c1 = 85.6917,
    c2 = -71567,
    w_ref = 477060,
    z = -1,
    a_0 = 3.72,
    MM = 0.12797) "HSeO3-";

  //303
  constant DataRecord HSeO4_n1(
    R = Modelica.Constants.R/HSeO4_n1.MM,
    G_ref = -452290,
    H_ref = -581576,
    S_ref = 149.369,
    a1 = 3.1482e-05,
    a2 = 4431.78,
    a3 = 6.6262e-05,
    a4 = -134591,
    c1 = 97.4265,
    c2 = -23836,
    w_ref = 455430,
    z = -1,
    a_0 = 3.72,
    MM = 0.14397) "HSeO4-";

  //304
  constant DataRecord HSiO3_n1(
    R = Modelica.Constants.R/HSiO3_n1.MM,
    G_ref = -1015879,
    H_ref = -1145880,
    S_ref = 20.92,
    a1 = 1.2441e-05,
    a2 = -215.81,
    a3 = 0.00024881,
    a4 = -115374,
    c1 = 34.095,
    c2 = -305947,
    w_ref = 648980,
    z = -1,
    a_0 = 3.72,
    MM = 0.077088) "HSiO3-";

  //305
  constant DataRecord HVO4_n2(
    R = Modelica.Constants.R/HVO4_n2.MM,
    G_ref = -974872,
    H_ref = -1158968,
    S_ref = 16.736,
    a1 = 1.8999e-05,
    a2 = 1383.94,
    a3 = 0.00018606,
    a4 = -121989,
    c1 = 361.3951,
    c2 = 617073,
    w_ref = 1319090,
    z = -2,
    a_0 = 3.72,
    MM = 0.11595) "HVO4-2";

  //306
  constant DataRecord He_aq(
    R = Modelica.Constants.R/He_aq.MM,
    G_ref = 19489,
    H_ref = -628,
    S_ref = 58.66,
    a1 = 1.4528e-05,
    a2 = 291.5,
    a3 = 0.00022912,
    a4 = -117478,
    c1 = 109.0798,
    c2 = 190933,
    w_ref = -88830,
    z = 0,
    a_0 = 3.72,
    MM = 0.004) "He,aq";

  //307
  constant DataRecord Hg_p2(
    R = Modelica.Constants.R/Hg_p2.MM,
    G_ref = 164682,
    H_ref = 170163,
    S_ref = -36.317,
    a1 = -2.2092e-06,
    a2 = -3795.18,
    a3 = 0.00038975,
    a4 = -100583,
    c1 = 75.5685,
    c2 = -108219,
    w_ref = 481660,
    z = 2,
    a_0 = 5,
    MM = 0.20059) "Hg+2";

  //308
  constant DataRecord HgOH_p1(
    R = Modelica.Constants.R/HgOH_p1.MM,
    G_ref = -53137,
    H_ref = -85354,
    S_ref = 69.873,
    a1 = 5.2321e-06,
    a2 = -1974.68,
    a3 = 0.00031745,
    a4 = -108106,
    c1 = 151.3805,
    c2 = 269341,
    w_ref = 125440,
    z = 1,
    a_0 = 3.72,
    MM = 0.2176) "HgOH+";

  //309
  constant DataRecord HgO_aq(
    R = Modelica.Constants.R/HgO_aq.MM,
    G_ref = -37238,
    H_ref = -80333,
    S_ref = 34.309,
    a1 = 2.5903e-06,
    a2 = -2621.86,
    a3 = 0.0003433,
    a4 = -105428,
    c1 = -36.1665,
    c2 = -338331,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.21659) "HgO,aq";

  //310
  constant DataRecord HHgO2_n1(
    R = Modelica.Constants.R/HHgO2_n1.MM,
    G_ref = -189117,
    H_ref = -308779,
    S_ref = 9.623,
    a1 = 7.5806e-06,
    a2 = -1404.65,
    a3 = 0.00029576,
    a4 = -110462,
    c1 = 15.7268,
    c2 = -375828,
    w_ref = 667850,
    z = -1,
    a_0 = 3.72,
    MM = 0.2336) "HHgO2-";

  //311
  constant DataRecord Hg2_p2(
    R = Modelica.Constants.R/Hg2_p2.MM,
    G_ref = 153595,
    H_ref = 166858,
    S_ref = 65.731,
    a1 = 1.6846e-05,
    a2 = 856.84,
    a3 = 0.00020707,
    a4 = -119813,
    c1 = 91.2489,
    c2 = -9347,
    w_ref = 343130,
    z = 2,
    a_0 = 4,
    MM = 0.40118) "Hg2+2";

  //312
  constant DataRecord Ho_p2(
    R = Modelica.Constants.R/Ho_p2.MM,
    G_ref = -405430,
    H_ref = -394133,
    S_ref = -17.991,
    a1 = 3.93e-08,
    a2 = -3243.27,
    a3 = 0.0003674,
    a4 = -102859,
    c1 = 41.8174,
    c2 = -221568,
    w_ref = 469280,
    z = 2,
    a_0 = 3.72,
    MM = 0.16493) "Ho+2";

  //313
  constant DataRecord Ho_p3(
    R = Modelica.Constants.R/Ho_p3.MM,
    G_ref = -675298,
    H_ref = -707096,
    S_ref = -227.191,
    a1 = -1.3053e-05,
    a2 = -6443.03,
    a3 = 0.00049382,
    a4 = -89638,
    c1 = 36.2694,
    c2 = -410777,
    w_ref = 999930,
    z = 3,
    a_0 = 3.72,
    MM = 0.16493) "Ho+3";

  //314
  constant DataRecord Ho_p4(
    R = Modelica.Constants.R/Ho_p4.MM,
    G_ref = -127612,
    H_ref = -202087,
    S_ref = -435.973,
    a1 = -1.806e-05,
    a2 = -7662.74,
    a3 = 0.00054112,
    a4 = -84592,
    c1 = 167.3445,
    c2 = -137189,
    w_ref = 1568330,
    z = 4,
    a_0 = 3.72,
    MM = 0.16493) "Ho+4";

  //315
  constant DataRecord I_n1(
    R = Modelica.Constants.R/I_n1.MM,
    G_ref = -51923,
    H_ref = -56902,
    S_ref = 106.692,
    a1 = 3.2477e-05,
    a2 = 3462.76,
    a3 = 6.1124e-05,
    a4 = -130587,
    c1 = -26.2337,
    c2 = -206857,
    w_ref = 541160,
    z = -1,
    a_0 = 3,
    MM = 0.1269) "I-";

  //316
  constant DataRecord I3_n1(
    R = Modelica.Constants.R/I3_n1.MM,
    G_ref = -51463,
    H_ref = -51463,
    S_ref = 239.325,
    a1 = 4.1212e-05,
    a2 = 6807.24,
    a3 = -2.6974e-05,
    a4 = -144415,
    c1 = 87.3251,
    c2 = -15318,
    w_ref = 319160,
    z = -1,
    a_0 = 3.72,
    MM = 0.3807) "I3-";

  //317
  constant DataRecord HIO_aq(
    R = Modelica.Constants.R/HIO_aq.MM,
    G_ref = -99161,
    H_ref = -138072,
    S_ref = 95.395,
    a1 = 1.6309e-05,
    a2 = 726.09,
    a3 = 0.0002121,
    a4 = -119269,
    c1 = 3.0376,
    c2 = -205372,
    w_ref = -2220,
    z = 0,
    a_0 = 3.72,
    MM = 0.14391) "HIO,aq";

  //318
  constant DataRecord IO_n1(
    R = Modelica.Constants.R/IO_n1.MM,
    G_ref = -38493,
    H_ref = -107529,
    S_ref = -5.858,
    a1 = 2.5652e-06,
    a2 = -2625.75,
    a3 = 0.00034302,
    a4 = -105416,
    c1 = 4.6744,
    c2 = -421852,
    w_ref = 691620,
    z = -1,
    a_0 = 3.72,
    MM = 0.1429) "IO-";

  //319
  constant DataRecord IO3_n1(
    R = Modelica.Constants.R/IO3_n1.MM,
    G_ref = -128030,
    H_ref = -221334,
    S_ref = 118.407,
    a1 = 2.3911e-05,
    a2 = 2582.57,
    a3 = 0.00013908,
    a4 = -126951,
    c1 = 32.3394,
    c2 = -265035,
    w_ref = 502160,
    z = -1,
    a_0 = 4.25,
    MM = 0.1749) "IO3-";

  //320
  constant DataRecord IO4_n1(
    R = Modelica.Constants.R/IO4_n1.MM,
    G_ref = -58576,
    H_ref = -151461,
    S_ref = 221.752,
    a1 = 3.884e-05,
    a2 = 6226.92,
    a3 = -3.904e-06,
    a4 = -142013,
    c1 = 52.7523,
    c2 = -144009,
    w_ref = 345770,
    z = -1,
    a_0 = 3.5,
    MM = 0.1909) "IO4-";

  //321
  constant DataRecord ISOLEUCINE_aq(
    R = Modelica.Constants.R/ISOLEUCINE_aq.MM,
    G_ref = -343925,
    H_ref = -631366,
    S_ref = 220.915,
    a1 = 6.7109e-05,
    a2 = 10017.08,
    a3 = 0.00051467,
    a4 = -157678,
    c1 = 347.359,
    c2 = 89462,
    w_ref = -192300,
    z = 0,
    a_0 = 3.72,
    MM = 0.13117) "ISOLEUCINE,aq";

  //322
  constant DataRecord In_p3(
    R = Modelica.Constants.R/In_p3.MM,
    G_ref = -97906,
    H_ref = -104600,
    S_ref = -263.592,
    a1 = -9.3558e-06,
    a2 = -5538.7,
    a3 = 0.00045787,
    a4 = -93370,
    c1 = 85.2465,
    c2 = -259918,
    w_ref = 1060480,
    z = 3,
    a_0 = 9,
    MM = 0.11482) "In+3";

  //323
  constant DataRecord InOH_p2(
    R = Modelica.Constants.R/InOH_p2.MM,
    G_ref = -312126,
    H_ref = -396225,
    S_ref = -187.025,
    a1 = 8.5236e-06,
    a2 = -1173.95,
    a3 = 0.00028658,
    a4 = -111416,
    c1 = 126.579,
    c2 = -9347,
    w_ref = 726590,
    z = 2,
    a_0 = 3.72,
    MM = 0.13183) "InOH+2";

  //324
  constant DataRecord InO_p1(
    R = Modelica.Constants.R/InO_p1.MM,
    G_ref = -290370,
    H_ref = -321750,
    S_ref = -10.878,
    a1 = 1.0796e-05,
    a2 = -620.61,
    a3 = 0.0002652,
    a4 = -113704,
    c1 = -93.1584,
    c2 = -619583,
    w_ref = 247150,
    z = 1,
    a_0 = 3.72,
    MM = 0.13082) "InO+";

  //325
  constant DataRecord HInO2_aq(
    R = Modelica.Constants.R/HInO2_aq.MM,
    G_ref = -501243,
    H_ref = -565258,
    S_ref = 113.805,
    a1 = 1.7934e-05,
    a2 = 1122.61,
    a3 = 0.00019659,
    a4 = -120909,
    c1 = -186.2344,
    c2 = -859925,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14783) "HInO2,aq";

  //326
  constant DataRecord InO2_n1(
    R = Modelica.Constants.R/InO2_n1.MM,
    G_ref = -446433,
    H_ref = -524255,
    S_ref = 66.526,
    a1 = 1.8339e-05,
    a2 = 1222.73,
    a3 = 0.00019239,
    a4 = -121323,
    c1 = -85.6402,
    c2 = -700548,
    w_ref = 581620,
    z = -1,
    a_0 = 3.72,
    MM = 0.14682) "InO2-";

  //327
  constant DataRecord K_p1(
    R = Modelica.Constants.R/K_p1.MM,
    G_ref = -282462,
    H_ref = -252170,
    S_ref = 101.044,
    a1 = 1.4891e-05,
    a2 = -616.3,
    a3 = 0.0002274,
    a4 = -113470,
    c1 = 30.9616,
    c2 = -74935,
    w_ref = 80630,
    z = 1,
    a_0 = 3,
    MM = 0.0391) "K+";

  //328
  constant DataRecord KBr_aq(
    R = Modelica.Constants.R/KBr_aq.MM,
    G_ref = -376602,
    H_ref = -361163,
    S_ref = 198.74,
    a1 = 3.3715e-05,
    a2 = 4977.7,
    a3 = 4.4643e-05,
    a4 = -136850,
    c1 = 6.0651,
    c2 = -194891,
    w_ref = -2090,
    z = 0,
    a_0 = 3.72,
    MM = 0.11901) "KBr,aq";

  //329
  constant DataRecord KCl_aq(
    R = Modelica.Constants.R/KCl_aq.MM,
    G_ref = -399279,
    H_ref = -399112,
    S_ref = 176.774,
    a1 = 2.9259e-05,
    a2 = 3889.86,
    a3 = 8.74e-05,
    a4 = -132352,
    c1 = 3.984,
    c2 = -197707,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.07455) "KCl,aq";

  //330
  constant DataRecord KAlO2_aq(
    R = Modelica.Constants.R/KAlO2_aq.MM,
    G_ref = -1106170,
    H_ref = -1150525,
    S_ref = 149.369,
    a1 = 2.5133e-05,
    a2 = 2881.02,
    a3 = 0.00012734,
    a4 = -128185,
    c1 = 5.8325,
    c2 = 50710,
    w_ref = -20920,
    z = 0,
    a_0 = 3.72,
    MM = 0.09808) "KAlO2,aq";

  //331
  constant DataRecord KHSO4_aq(
    R = Modelica.Constants.R/KHSO4_aq.MM,
    G_ref = -1018386,
    H_ref = -1118802,
    S_ref = 230.12,
    a1 = 3.8169e-05,
    a2 = 6065.29,
    a3 = 1.895e-06,
    a4 = -141344,
    c1 = 167.2968,
    c2 = 364970,
    w_ref = -420,
    z = 0,
    a_0 = 3.72,
    MM = 0.13617) "KHSO4,aq";

  //332
  constant DataRecord KI_aq(
    R = Modelica.Constants.R/KI_aq.MM,
    G_ref = -325264,
    H_ref = -300533,
    S_ref = 205.853,
    a1 = 4.1158e-05,
    a2 = 6795.07,
    a3 = -2.6786e-05,
    a4 = -144361,
    c1 = 11.3863,
    c2 = -176397,
    w_ref = -2090,
    z = 0,
    a_0 = 3.72,
    MM = 0.166) "KI,aq";

  //333
  constant DataRecord KOH_aq(
    R = Modelica.Constants.R/KOH_aq.MM,
    G_ref = -437228,
    H_ref = -502917,
    S_ref = 108.366,
    a1 = 1.5873e-05,
    a2 = 620.86,
    a3 = 0.00021597,
    a4 = -118834,
    c1 = -25.6228,
    c2 = -301683,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.056108) "KOH,aq";

  //334
  constant DataRecord KSO4_n1(
    R = Modelica.Constants.R/KSO4_n1.MM,
    G_ref = -1031942,
    H_ref = -1158884,
    S_ref = 146.44,
    a1 = 2.4856e-05,
    a2 = 2814.74,
    a3 = 0.00012966,
    a4 = -127909,
    c1 = 41.4588,
    c2 = -219865,
    w_ref = 460070,
    z = -1,
    a_0 = 3.72,
    MM = 0.13516) "KSO4-";

  //335
  constant DataRecord Kr_aq(
    R = Modelica.Constants.R/Kr_aq.MM,
    G_ref = 14870,
    H_ref = -15272,
    S_ref = 63.011,
    a1 = 2.6242e-05,
    a2 = 3151.93,
    a3 = 0.0001167,
    a4 = -129302,
    c1 = 158.2498,
    c2 = 363945,
    w_ref = -95440,
    z = 0,
    a_0 = 3.72,
    MM = 0.0838) "Kr,aq";

  //336
  constant DataRecord LACTATE_aq(
    R = Modelica.Constants.R/LACTATE_aq.MM,
    G_ref = -512666,
    H_ref = -686469,
    S_ref = 133.47,
    a1 = 4.1212e-05,
    a2 = 5365.73,
    a3 = 0.00033886,
    a4 = -138453,
    c1 = 188.0783,
    c2 = -158251,
    w_ref = 479860,
    z = -1,
    a_0 = 3.72,
    MM = 0.089069) "LACTATE,aq";

  //337
  constant DataRecord LACTIC_ACID_aq(
    R = Modelica.Constants.R/LACTIC_ACID_aq.MM,
    G_ref = -534715,
    H_ref = -686176,
    S_ref = 208.363,
    a1 = 4.6744e-05,
    a2 = 7214.55,
    a3 = 0.00015934,
    a4 = -146093,
    c1 = 255.3073,
    c2 = 68157,
    w_ref = -107610,
    z = 0,
    a_0 = 3.72,
    MM = 0.090077) "LACTIC-ACID,aq";

  //338
  constant DataRecord LEUCINE_aq(
    R = Modelica.Constants.R/LEUCINE_aq.MM,
    G_ref = -343088,
    H_ref = -632077,
    S_ref = 215.476,
    a1 = 6.8351e-05,
    a2 = 10263.85,
    a3 = 0.00051709,
    a4 = -158699,
    c1 = 345.4273,
    c2 = 174766,
    w_ref = -184050,
    z = 0,
    a_0 = 3.72,
    MM = 0.13117) "LEUCINE,aq";

  //339
  constant DataRecord La_p3(
    R = Modelica.Constants.R/La_p3.MM,
    G_ref = -686176,
    H_ref = -709606,
    S_ref = -217.568,
    a1 = -1.1665e-05,
    a2 = -6017.6,
    a3 = 0.00045858,
    a4 = -91395,
    c1 = 17.7376,
    c2 = -444014,
    w_ref = 902570,
    z = 3,
    a_0 = 9,
    MM = 0.13891) "La+3";

  //340
  constant DataRecord Li_p1(
    R = Modelica.Constants.R/Li_p1.MM,
    G_ref = -292600,
    H_ref = -278454,
    S_ref = 11.297,
    a1 = -9.92e-08,
    a2 = -28.87,
    a3 = 0.00048451,
    a4 = -116152,
    c1 = 80.3328,
    c2 = -10042,
    w_ref = 203430,
    z = 1,
    a_0 = 6,
    MM = 0.00694) "Li+";

  //341
  constant DataRecord LiCl_aq(
    R = Modelica.Constants.R/LiCl_aq.MM,
    G_ref = -415262,
    H_ref = -442165,
    S_ref = 54.978,
    a1 = 2.3362e-05,
    a2 = 2449.9,
    a3 = 0.000144,
    a4 = -126399,
    c1 = 74.1137,
    c2 = 46049,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.04239) "LiCl,aq";

  //342
  constant DataRecord LiOH_aq(
    R = Modelica.Constants.R/LiOH_aq.MM,
    G_ref = -451872,
    H_ref = -508356,
    S_ref = 7.95,
    a1 = 9.5182e-06,
    a2 = -930.61,
    a3 = 0.0002769,
    a4 = -112424,
    c1 = 41.3178,
    c2 = -69007,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.023948) "LiOH,aq";

  //343
  constant DataRecord Lu_p3(
    R = Modelica.Constants.R/Lu_p3.MM,
    G_ref = -666930,
    H_ref = -702494,
    S_ref = -264.01,
    a1 = -1.4908e-05,
    a2 = -6895.73,
    a3 = 0.00051161,
    a4 = -87768,
    c1 = 40.02,
    c2 = -406517,
    w_ref = 1027340,
    z = 3,
    a_0 = 3.72,
    MM = 0.17496) "Lu+3";

  //344
  constant DataRecord Lu_p4(
    R = Modelica.Constants.R/Lu_p4.MM,
    G_ref = 112968,
    H_ref = 41003,
    S_ref = -451.872,
    a1 = -1.8195e-05,
    a2 = -7695.25,
    a3 = 0.00054228,
    a4 = -84458,
    c1 = 166.9592,
    c2 = -147419,
    w_ref = 1596070,
    z = 4,
    a_0 = 3.72,
    MM = 0.17496) "Lu+4";

  //345
  constant DataRecord m_TOLUATE_aq(
    R = Modelica.Constants.R/m_TOLUATE_aq.MM,
    G_ref = -218196,
    H_ref = -398944,
    S_ref = 166.942,
    a1 = 6.6155e-05,
    a2 = 12179.58,
    a3 = -8.4098e-05,
    a4 = -166619,
    c1 = 402.9037,
    c2 = -130968,
    w_ref = 428530,
    z = -1,
    a_0 = 3.72,
    MM = 0.13514) "m-TOLUATE,aq";

  //346
  constant DataRecord m_TOLUIC_ACID_aq(
    R = Modelica.Constants.R/m_TOLUIC_ACID_aq.MM,
    G_ref = -242505,
    H_ref = -399363,
    S_ref = 248.111,
    a1 = 7.112e-05,
    a2 = 13304.62,
    a3 = -0.00010964,
    a4 = -171272,
    c1 = 442.6572,
    c2 = 296754,
    w_ref = -81250,
    z = 0,
    a_0 = 3.72,
    MM = 0.13614) "m-TOLUIC-ACID,aq";

  //347
  constant DataRecord MALONATE_aq(
    R = Modelica.Constants.R/MALONATE_aq.MM,
    G_ref = -685172,
    H_ref = -874456,
    S_ref = 53.555,
    a1 = 3.3657e-05,
    a2 = 4812.86,
    a3 = 8.3442e-05,
    a4 = -136164,
    c1 = -115.5997,
    c2 = -204953,
    w_ref = 1263230,
    z = -2,
    a_0 = 3.72,
    MM = 0.10205) "MALONATE,aq";

  //348
  constant DataRecord MALONIC_ACID_aq(
    R = Modelica.Constants.R/MALONIC_ACID_aq.MM,
    G_ref = -733957,
    H_ref = -869728,
    S_ref = 233.049,
    a1 = 4.5333e-05,
    a2 = 7461.75,
    a3 = 2.2744e-05,
    a4 = -147118,
    c1 = 162.1802,
    c2 = -48827,
    w_ref = -91250,
    z = 0,
    a_0 = 3.72,
    MM = 0.10406) "MALONIC-ACID,aq";

  //349
  constant DataRecord METHANAMINE_aq(
    R = Modelica.Constants.R/METHANAMINE_aq.MM,
    G_ref = 21087,
    H_ref = -68283,
    S_ref = 127.612,
    a1 = 3.119e-05,
    a2 = 2837.76,
    a3 = 0.00045557,
    a4 = -128001,
    c1 = 163.3116,
    c2 = -64965,
    w_ref = -51000,
    z = 0,
    a_0 = 3.72,
    MM = 0.031056) "METHANAMINE,aq";

  //350
  constant DataRecord METHANE_aq(
    R = Modelica.Constants.R/METHANE_aq.MM,
    G_ref = -34451,
    H_ref = -87906,
    S_ref = 87.822,
    a1 = 2.8291e-05,
    a2 = 3651.75,
    a3 = 9.7119e-05,
    a4 = -131365,
    c1 = 176.1217,
    c2 = 438094,
    w_ref = -133010,
    z = 0,
    a_0 = 3.72,
    MM = 0.016042) "METHANE,aq";

  //351
  constant DataRecord METHANOL_aq(
    R = Modelica.Constants.R/METHANOL_aq.MM,
    G_ref = -175937,
    H_ref = -246312,
    S_ref = 134.725,
    a1 = 2.903e-05,
    a2 = 2307.31,
    a3 = 0.00047705,
    a4 = -125809,
    c1 = 165.2061,
    c2 = -62701,
    w_ref = -61760,
    z = 0,
    a_0 = 3.72,
    MM = 0.032042) "METHANOL,aq";

  //352
  constant DataRecord METHIONINE_aq(
    R = Modelica.Constants.R/METHIONINE_aq.MM,
    G_ref = -502917,
    H_ref = -743078,
    S_ref = 274.889,
    a1 = 6.6747e-05,
    a2 = 9944.87,
    a3 = 0.00051406,
    a4 = -157381,
    c1 = 261.8694,
    c2 = 28368,
    w_ref = -274050,
    z = 0,
    a_0 = 3.72,
    MM = 0.1492) "METHIONINE,aq";

  //353
  constant DataRecord Mg_CO3_aq(
    R = Modelica.Constants.R/Mg_CO3_aq.MM,
    G_ref = -998972,
    H_ref = -1132065,
    S_ref = -100.416,
    a1 = -2.2803e-06,
    a2 = -3812.88,
    a3 = 0.00039045,
    a4 = -100500,
    c1 = -43.9278,
    c2 = -364259,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.08432) "Mg(CO3),aq";

  //354
  constant DataRecord Mg_HCO3_p1(
    R = Modelica.Constants.R/Mg_HCO3_p1.MM,
    G_ref = -1046837,
    H_ref = -1153738,
    S_ref = -12.552,
    a1 = 1.3686e-05,
    a2 = 86.19,
    a3 = 0.00023719,
    a4 = -116650,
    c1 = 197.8363,
    c2 = 390786,
    w_ref = 250620,
    z = 1,
    a_0 = 3.72,
    MM = 0.085328) "Mg(HCO3)+";

  //355
  constant DataRecord Mg_p2(
    R = Modelica.Constants.R/Mg_p2.MM,
    G_ref = -453985,
    H_ref = -465960,
    S_ref = -138.072,
    a1 = -3.438e-06,
    a2 = -3597.82,
    a3 = 0.00035104,
    a4 = -99998,
    c1 = 87.0272,
    c2 = -246521,
    w_ref = 643160,
    z = 2,
    a_0 = 8,
    MM = 0.02431) "Mg+2";

  //356
  constant DataRecord MgCl_p1(
    R = Modelica.Constants.R/MgCl_p1.MM,
    G_ref = -584505,
    H_ref = -631504,
    S_ref = -79.496,
    a1 = 9.301e-06,
    a2 = -983.45,
    a3 = 0.00027894,
    a4 = -112207,
    c1 = 119.6691,
    c2 = 86107,
    w_ref = 353510,
    z = 1,
    a_0 = 3.72,
    MM = 0.05976) "MgCl+";

  //357
  constant DataRecord MgF_p1(
    R = Modelica.Constants.R/MgF_p1.MM,
    G_ref = -743455,
    H_ref = -798935,
    S_ref = -117.445,
    a1 = -1.2447e-06,
    a2 = -3558.45,
    a3 = 0.00038015,
    a4 = -101562,
    c1 = 159.092,
    c2 = 206275,
    w_ref = 406100,
    z = 1,
    a_0 = 3.72,
    MM = 0.04331) "MgF+";

  //358
  constant DataRecord MgOH_p1(
    R = Modelica.Constants.R/MgOH_p1.MM,
    G_ref = -624483,
    H_ref = -639964,
    S_ref = -79.914,
    a1 = 9.6671e-06,
    a2 = -893.91,
    a3 = 0.00027542,
    a4 = -112575,
    c1 = 133.8913,
    c2 = 135536,
    w_ref = 353510,
    z = 1,
    a_0 = 3.72,
    MM = 0.041318) "MgOH+";

  //359
  constant DataRecord Mn_p2(
    R = Modelica.Constants.R/Mn_p2.MM,
    G_ref = -230538,
    H_ref = -221334,
    S_ref = -67.781,
    a1 = -4.247e-07,
    a2 = -3348.87,
    a3 = 0.00036987,
    a4 = -102424,
    c1 = 69.7356,
    c2 = -161908,
    w_ref = 586010,
    z = 2,
    a_0 = 6,
    MM = 0.05494) "Mn+2";

  //360
  constant DataRecord Mn_p3(
    R = Modelica.Constants.R/Mn_p3.MM,
    G_ref = -84935,
    H_ref = -128449,
    S_ref = -309.616,
    a1 = -1.2267e-05,
    a2 = -6248.39,
    a3 = 0.00048552,
    a4 = -90437,
    c1 = 67.1996,
    c2 = -345147,
    w_ref = 1130730,
    z = 3,
    a_0 = 3.72,
    MM = 0.05494) "Mn+3";

  //361
  constant DataRecord MnCl_p1(
    R = Modelica.Constants.R/MnCl_p1.MM,
    G_ref = -361037,
    H_ref = -369364,
    S_ref = 50.208,
    a1 = 1.1765e-05,
    a2 = -381.83,
    a3 = 0.0002553,
    a4 = -114692,
    c1 = 101.6034,
    c2 = 87128,
    w_ref = 154220,
    z = 1,
    a_0 = 3.72,
    MM = 0.09039) "MnCl+";

  //362
  constant DataRecord MnOH_p1(
    R = Modelica.Constants.R/MnOH_p1.MM,
    G_ref = -407103,
    H_ref = -446851,
    S_ref = 1.255,
    a1 = 1.3443e-06,
    a2 = -2924.66,
    a3 = 0.00035489,
    a4 = -104177,
    c1 = 68.1967,
    c2 = -52815,
    w_ref = 228610,
    z = 1,
    a_0 = 3.72,
    MM = 0.071948) "MnOH+";

  //363
  constant DataRecord MnO_aq(
    R = Modelica.Constants.R/MnO_aq.MM,
    G_ref = -340996,
    H_ref = -414634,
    S_ref = -10.46,
    a1 = -1.573e-07,
    a2 = -3290.97,
    a3 = 0.00036915,
    a4 = -102663,
    c1 = -6.4969,
    c2 = -235204,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.07094) "MnO,aq";

  //364
  constant DataRecord HMnO2_n1(
    R = Modelica.Constants.R/HMnO2_n1.MM,
    G_ref = -506264,
    H_ref = -627182,
    S_ref = -38.074,
    a1 = 4.3321e-06,
    a2 = -2195.76,
    a3 = 0.00032637,
    a4 = -107190,
    c1 = 85.8787,
    c2 = -155088,
    w_ref = 739940,
    z = -1,
    a_0 = 3.72,
    MM = 0.087948) "HMnO2-";

  //365
  constant DataRecord MnO2_n2(
    R = Modelica.Constants.R/MnO2_n2.MM,
    G_ref = -429278,
    H_ref = -557727,
    S_ref = -63.597,
    a1 = 4.807e-06,
    a2 = -2079.99,
    a3 = 0.0003219,
    a4 = -107671,
    c1 = -16.8929,
    c2 = -736346,
    w_ref = 1439630,
    z = -2,
    a_0 = 3.72,
    MM = 0.08694) "MnO2-2";

  //366
  constant DataRecord MnO4_n1(
    R = Modelica.Constants.R/MnO4_n1.MM,
    G_ref = -450198,
    H_ref = -543502,
    S_ref = 194.556,
    a1 = 3.2739e-05,
    a2 = 4739.51,
    a3 = 5.4024e-05,
    a4 = -135863,
    c1 = 57.035,
    c2 = -142306,
    w_ref = 386940,
    z = -1,
    a_0 = 3.5,
    MM = 0.11894) "MnO4-";

  //367
  constant DataRecord MnO4_n2(
    R = Modelica.Constants.R/MnO4_n2.MM,
    G_ref = -503754,
    H_ref = -655214,
    S_ref = 64.852,
    a1 = 2.368e-05,
    a2 = 2525.8,
    a3 = 0.00014136,
    a4 = -126712,
    c1 = -22.1375,
    c2 = -692879,
    w_ref = 1246960,
    z = -2,
    a_0 = 3.72,
    MM = 0.11894) "MnO4-2";

  //368
  constant DataRecord MnSO4_aq(
    R = Modelica.Constants.R/MnSO4_aq.MM,
    G_ref = -985918,
    H_ref = -1116082,
    S_ref = 20.92,
    a1 = 8.9345e-06,
    a2 = -1072.99,
    a3 = 0.00028246,
    a4 = -111834,
    c1 = -26.1768,
    c2 = -302537,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.151) "MnSO4,aq";

  //369
  constant DataRecord HMoO4_n1(
    R = Modelica.Constants.R/HMoO4_n1.MM,
    G_ref = -863578,
    H_ref = -997884,
    S_ref = 135.98,
    a1 = 3.2351e-05,
    a2 = 4642.82,
    a3 = 5.8191e-05,
    a4 = -135461,
    c1 = 139.7163,
    c2 = 116784,
    w_ref = 475260,
    z = -1,
    a_0 = 3.72,
    MM = 0.16095) "HMoO4-";

  //370
  constant DataRecord MoO4_n2(
    R = Modelica.Constants.R/MoO4_n2.MM,
    G_ref = -838474,
    H_ref = -997047,
    S_ref = 37.656,
    a1 = 2.9142e-05,
    a2 = 1133.65,
    a3 = 0.00078081,
    a4 = -120955,
    c1 = 27.9613,
    c2 = -531795,
    w_ref = 1287710,
    z = -2,
    a_0 = 4.5,
    MM = 0.15994) "MoO4-2";

  //371
  constant DataRecord N_BUTANE_aq(
    R = Modelica.Constants.R/N_BUTANE_aq.MM,
    G_ref = 151,
    H_ref = -151586,
    S_ref = 167.444,
    a1 = 5.3934e-05,
    a2 = 9914.41,
    a3 = -0.0001493,
    a4 = -157256,
    c1 = 330.7741,
    c2 = 1014235,
    w_ref = -253590,
    z = 0,
    a_0 = 3.72,
    MM = 0.058119) "N-BUTANE,aq";

  //372
  constant DataRecord BUTANOATE_aq(
    R = Modelica.Constants.R/BUTANOATE_aq.MM,
    G_ref = -354176,
    H_ref = -538188,
    S_ref = 133.051,
    a1 = 4.9256e-05,
    a2 = 7133.39,
    a3 = 0.00031153,
    a4 = -145758,
    c1 = 260.7197,
    c2 = -149227,
    w_ref = 479860,
    z = -1,
    a_0 = 3.72,
    MM = 0.087095) "BUTANOATE,aq";

  //373
  constant DataRecord N_BUTYLBENZENE_aq(
    R = Modelica.Constants.R/N_BUTYLBENZENE_aq.MM,
    G_ref = 151084,
    H_ref = -57823,
    S_ref = 263.592,
    a1 = 9.4113e-05,
    a2 = 18516.42,
    a3 = -0.00022812,
    a4 = -192815,
    c1 = 555.2917,
    c2 = 434822,
    w_ref = -71000,
    z = 0,
    a_0 = 3.72,
    MM = 0.13421) "N-BUTYLBENZENE,aq";

  //374
  constant DataRecord N_HEPTANE_aq(
    R = Modelica.Constants.R/N_HEPTANE_aq.MM,
    G_ref = 27070,
    H_ref = -221543,
    S_ref = 251.04,
    a1 = 8.0644e-05,
    a2 = 16435.17,
    a3 = -0.00040534,
    a4 = -184213,
    c1 = 472.3652,
    c2 = 1546913,
    w_ref = -380160,
    z = 0,
    a_0 = 3.72,
    MM = 0.1002) "N-HEPTANE,aq";

  //375
  constant DataRecord N_HEPTYLBENZENE_aq(
    R = Modelica.Constants.R/N_HEPTYLBENZENE_aq.MM,
    G_ref = 181376,
    H_ref = -516389,
    S_ref = 347.69,
    a1 = 0.00012087,
    a2 = 24580.96,
    a3 = -0.00036609,
    a4 = -217886,
    c1 = 731.9736,
    c2 = 646888,
    w_ref = -15230,
    z = 0,
    a_0 = 3.72,
    MM = 0.17629) "N-HEPTYLBENZENE,aq";

  //376
  constant DataRecord N_HEXANE_aq(
    R = Modelica.Constants.R/N_HEXANE_aq.MM,
    G_ref = 18493,
    H_ref = -198322,
    S_ref = 221.334,
    a1 = 7.175e-05,
    a2 = 14264.39,
    a3 = -0.0003202,
    a4 = -175238,
    c1 = 424.5262,
    c2 = 1366227,
    w_ref = -335180,
    z = 0,
    a_0 = 3.72,
    MM = 0.086171) "N-HEXANE,aq";

  //377
  constant DataRecord N_HEXYLBENZENE_aq(
    R = Modelica.Constants.R/N_HEXYLBENZENE_aq.MM,
    G_ref = 171042,
    H_ref = -104558,
    S_ref = 319.658,
    a1 = 0.00011142,
    a2 = 22438.29,
    a3 = -0.00031734,
    a4 = -209028,
    c1 = 673.0797,
    c2 = 576200,
    w_ref = -33850,
    z = 0,
    a_0 = 3.72,
    MM = 0.16226) "N-HEXYLBENZENE,aq";

  //378
  constant DataRecord N_OCTANE_aq(
    R = Modelica.Constants.R/N_OCTANE_aq.MM,
    G_ref = 35899,
    H_ref = -248571,
    S_ref = 266.939,
    a1 = 8.9609e-05,
    a2 = 18624.07,
    a3 = -0.00049136,
    a4 = -193263,
    c1 = 522.1314,
    c2 = 1727595,
    w_ref = -404260,
    z = 0,
    a_0 = 3.72,
    MM = 0.11422) "N-OCTANE,aq";

  //379
  constant DataRecord N_OCTYLBENZENE_aq(
    R = Modelica.Constants.R/N_OCTYLBENZENE_aq.MM,
    G_ref = 189954,
    H_ref = -579066,
    S_ref = 375.723,
    a1 = 0.00013032,
    a2 = 26723.71,
    a3 = -0.00041484,
    a4 = -226744,
    c1 = 790.868,
    c2 = 717573,
    w_ref = 3350,
    z = 0,
    a_0 = 3.72,
    MM = 0.19031) "N-OCTYLBENZENE,aq";

  //380
  constant DataRecord N_PENTANE_aq(
    R = Modelica.Constants.R/N_PENTANE_aq.MM,
    G_ref = 8912,
    H_ref = -173887,
    S_ref = 198.74,
    a1 = 6.282e-05,
    a2 = 12082.18,
    a3 = -0.00023409,
    a4 = -166218,
    c1 = 373.2421,
    c2 = 1177022,
    w_ref = -300960,
    z = 0,
    a_0 = 3.72,
    MM = 0.072145) "N-PENTANE,aq";

  //381
  constant DataRecord N_PENTYLBENZENE_aq(
    R = Modelica.Constants.R/N_PENTYLBENZENE_aq.MM,
    G_ref = 162297,
    H_ref = -81170,
    S_ref = 291.625,
    a1 = 0.00010196,
    a2 = 20295.66,
    a3 = -0.0002686,
    a4 = -200171,
    c1 = 614.1857,
    c2 = 505511,
    w_ref = -52430,
    z = 0,
    a_0 = 3.72,
    MM = 0.14824) "N-PENTYLBENZENE,aq";

  //382
  constant DataRecord N_PROPYLBENZENE_aq(
    R = Modelica.Constants.R/N_PROPYLBENZENE_aq.MM,
    G_ref = 143804,
    H_ref = -36108,
    S_ref = 231.375,
    a1 = 8.448e-05,
    a2 = 16332.33,
    a3 = -0.00017837,
    a4 = -183786,
    c1 = 517.8779,
    c2 = 391003,
    w_ref = -87030,
    z = 0,
    a_0 = 3.72,
    MM = 0.12018) "N-PROPYLBENZENE,aq";

  //383
  constant DataRecord N2_aq(
    R = Modelica.Constants.R/N2_aq.MM,
    G_ref = 18188,
    H_ref = -10439,
    S_ref = 95.814,
    a1 = 2.596e-05,
    a2 = 3082.98,
    a3 = 0.00011941,
    a4 = -129018,
    c1 = 149.75,
    c2 = 350310,
    w_ref = -145100,
    z = 0,
    a_0 = 3.72,
    MM = 0.028013) "N2,aq";

  //384
  constant DataRecord NH3_aq(
    R = Modelica.Constants.R/NH3_aq.MM,
    G_ref = -26706,
    H_ref = -81337,
    S_ref = 107.822,
    a1 = 2.1301e-05,
    a2 = 1170.26,
    a3 = 0.00036086,
    a4 = -121110,
    c1 = 84.9352,
    c2 = -48953,
    w_ref = -20920,
    z = 0,
    a_0 = 3.72,
    MM = 0.01703) "NH3,aq";

  //385
  constant DataRecord NH4_p1(
    R = Modelica.Constants.R/NH4_p1.MM,
    G_ref = -79454,
    H_ref = -133260,
    S_ref = 111.169,
    a1 = 1.6218e-05,
    a2 = 981.06,
    a3 = 0.00035817,
    a4 = -120328,
    c1 = 73.0108,
    c2 = -879,
    w_ref = 62840,
    z = 1,
    a_0 = 2.5,
    MM = 0.018038) "NH4+";

  //386
  constant DataRecord NO2_n1(
    R = Modelica.Constants.R/NO2_n1.MM,
    G_ref = -32217,
    H_ref = -104600,
    S_ref = 123.01,
    a1 = 2.3373e-05,
    a2 = 2451.41,
    a3 = 0.00014423,
    a4 = -126407,
    c1 = 14.3344,
    c2 = -325549,
    w_ref = 495680,
    z = -1,
    a_0 = 3,
    MM = 0.046007) "NO2-";

  //387
  constant DataRecord NO3_n1(
    R = Modelica.Constants.R/NO3_n1.MM,
    G_ref = -110905,
    H_ref = -206811,
    S_ref = 146.942,
    a1 = 3.0611e-05,
    a2 = 2837.76,
    a3 = -0.00019597,
    a4 = -128005,
    c1 = 32.2168,
    c2 = -281374,
    w_ref = 459280,
    z = -1,
    a_0 = 3,
    MM = 0.062007) "NO3-";

  //388
  constant DataRecord NONANOATE_aq(
    R = Modelica.Constants.R/NONANOATE_aq.MM,
    G_ref = -312838,
    H_ref = -656846,
    S_ref = 273.634,
    a1 = 9.3711e-05,
    a2 = 18426,
    a3 = -0.00022628,
    a4 = -192443,
    c1 = 686.3622,
    c2 = -93902,
    w_ref = 267230,
    z = -1,
    a_0 = 3.72,
    MM = 0.15722) "NONANOATE,aq";

  //389
  constant DataRecord NONANOIC_ACID_aq(
    R = Modelica.Constants.R/NONANOIC_ACID_aq.MM,
    G_ref = -339824,
    H_ref = -654922,
    S_ref = 374.886,
    a1 = 0.00010083,
    a2 = 20038.26,
    a3 = -0.00026245,
    a4 = -199108,
    c1 = 673.1085,
    c2 = 572066,
    w_ref = 2760,
    z = 0,
    a_0 = 3.72,
    MM = 0.15823) "NONANOIC-ACID,aq";

  //390
  constant DataRecord Na_p1(
    R = Modelica.Constants.R/Na_p1.MM,
    G_ref = -261881,
    H_ref = -240300,
    S_ref = 58.409,
    a1 = 7.6944e-06,
    a2 = -956.04,
    a3 = 0.00013623,
    a4 = -114056,
    c1 = 76.0651,
    c2 = -124725,
    w_ref = 138320,
    z = 1,
    a_0 = 4.25,
    MM = 0.02299) "Na+";

  //391
  constant DataRecord NaAl_OH4_aq(
    R = Modelica.Constants.R/NaAl_OH4_aq.MM,
    G_ref = -1567372,
    H_ref = -1730435,
    S_ref = 204.179,
    a1 = 3.8186e-05,
    a2 = 6000.32,
    a3 = 4.69e-06,
    a4 = -141080,
    c1 = 254.0345,
    c2 = -587948,
    w_ref = 0,
    z = 0,
    a_0 = 3.72,
    MM = 0.118) "NaAl(OH)4,aq";

  //392
  constant DataRecord NaAlO2_aq(
    R = Modelica.Constants.R/NaAlO2_aq.MM,
    G_ref = -1088999,
    H_ref = -1160052,
    S_ref = 46.442,
    a1 = 1.9933e-05,
    a2 = 202.92,
    a3 = 0.00053007,
    a4 = -204911,
    c1 = 135.6871,
    c2 = 158992,
    w_ref = -41840,
    z = 0,
    a_0 = 3.72,
    MM = 0.08197) "NaAlO2,aq";

  //393
  constant DataRecord NaBr_aq(
    R = Modelica.Constants.R/NaBr_aq.MM,
    G_ref = -358192,
    H_ref = -354929,
    S_ref = 142.674,
    a1 = 2.563e-05,
    a2 = 3003.65,
    a3 = 0.00012223,
    a4 = -128687,
    c1 = 42.6458,
    c2 = -59040,
    w_ref = -29290,
    z = 0,
    a_0 = 3.72,
    MM = 0.1029) "NaBr,aq";

  //394
  constant DataRecord NaCl_aq(
    R = Modelica.Constants.R/NaCl_aq.MM,
    G_ref = -388735,
    H_ref = -402333,
    S_ref = 117.152,
    a1 = 2.1072e-05,
    a2 = 1890.71,
    a3 = 0.00016598,
    a4 = -124089,
    c1 = 45.1788,
    c2 = -54522,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.05844) "NaCl,aq";

  //395
  constant DataRecord NaF_aq(
    R = Modelica.Constants.R/NaF_aq.MM,
    G_ref = -537937,
    H_ref = -568438,
    S_ref = 50.208,
    a1 = 1.0601e-05,
    a2 = -666.18,
    a3 = 0.00026648,
    a4 = -113516,
    c1 = 51.7996,
    c2 = -31510,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.04199) "NaF,aq";

  //396
  constant DataRecord NaHSiO3_aq(
    R = Modelica.Constants.R/NaHSiO3_aq.MM,
    G_ref = -1288212,
    H_ref = -1397012,
    S_ref = 41.84,
    a1 = 1.4614e-05,
    a2 = 313.8,
    a3 = 0.00022796,
    a4 = -117570,
    c1 = 84.6821,
    c2 = 82780,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.10008) "NaHSiO3,aq";

  //397
  constant DataRecord NaI_aq(
    R = Modelica.Constants.R/NaI_aq.MM,
    G_ref = -305014,
    H_ref = -289119,
    S_ref = 157.737,
    a1 = 3.2003e-05,
    a2 = 4559.68,
    a3 = 6.1074e-05,
    a4 = -135122,
    c1 = 50.6268,
    c2 = -40543,
    w_ref = -420,
    z = 0,
    a_0 = 3.72,
    MM = 0.14989) "NaI,aq";

  //398
  constant DataRecord NaOH_aq(
    R = Modelica.Constants.R/NaOH_aq.MM,
    G_ref = -417982,
    H_ref = -469863,
    S_ref = 44.769,
    a1 = 9.3462e-06,
    a2 = -974.33,
    a3 = 0.000279,
    a4 = -112240,
    c1 = 16.7971,
    c2 = -154235,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.039998) "NaOH,aq";

  //399
  constant DataRecord Nd_p2(
    R = Modelica.Constants.R/Nd_p2.MM,
    G_ref = -419237,
    H_ref = -402082,
    S_ref = -2.092,
    a1 = 3.598e-07,
    a2 = -3165.36,
    a3 = 0.00036441,
    a4 = -103182,
    c1 = 37.9137,
    c2 = -227534,
    w_ref = 445550,
    z = 2,
    a_0 = 3.72,
    MM = 0.14424) "Nd+2";

  //400
  constant DataRecord Nd_p3(
    R = Modelica.Constants.R/Nd_p3.MM,
    G_ref = -671950,
    H_ref = -696636,
    S_ref = -207.108,
    a1 = -1.4103e-05,
    a2 = -6085.71,
    a3 = 0.00034816,
    a4 = -91115,
    c1 = 6.7931,
    c2 = -495151,
    w_ref = 943490,
    z = 3,
    a_0 = 9,
    MM = 0.14424) "Nd+3";

  //401
  constant DataRecord Nd_p4(
    R = Modelica.Constants.R/Nd_p4.MM,
    G_ref = -196648,
    H_ref = -263174,
    S_ref = -412.961,
    a1 = -1.7883e-05,
    a2 = -7620.24,
    a3 = 0.00053956,
    a4 = -84768,
    c1 = 168.5198,
    c2 = -122704,
    w_ref = 1535820,
    z = 4,
    a_0 = 3.72,
    MM = 0.14424) "Nd+4";

  //402
  constant DataRecord Ne_aq(
    R = Modelica.Constants.R/Ne_aq.MM,
    G_ref = 19100,
    H_ref = -3640,
    S_ref = 70.04,
    a1 = 1.8706e-05,
    a2 = 1311.77,
    a3 = 0.00018902,
    a4 = -121696,
    c1 = 113.3768,
    c2 = 211388,
    w_ref = -106060,
    z = 0,
    a_0 = 3.72,
    MM = 0.02018) "Ne,aq";

  //403
  constant DataRecord Ni_p2(
    R = Modelica.Constants.R/Ni_p2.MM,
    G_ref = -45606,
    H_ref = -53974,
    S_ref = -128.867,
    a1 = -7.0885e-06,
    a2 = -4986.53,
    a3 = 0.00043658,
    a4 = -95659,
    c1 = 55.1891,
    c2 = -226685,
    w_ref = 630400,
    z = 2,
    a_0 = 3.72,
    MM = 0.05871) "Ni+2";

  //404
  constant DataRecord NiCl_p1(
    R = Modelica.Constants.R/NiCl_p1.MM,
    G_ref = -171209,
    H_ref = -215058,
    S_ref = -71.128,
    a1 = 4.7359e-06,
    a2 = -2098.15,
    a3 = 0.00032275,
    a4 = -107596,
    c1 = 76.9283,
    c2 = -57932,
    w_ref = 339360,
    z = 1,
    a_0 = 3.72,
    MM = 0.09416) "NiCl+";

  //405
  constant DataRecord NiOH_p1(
    R = Modelica.Constants.R/NiOH_p1.MM,
    G_ref = -221124,
    H_ref = -283048,
    S_ref = -75.312,
    a1 = -2.9602e-06,
    a2 = -3976.93,
    a3 = 0.00039648,
    a4 = -99830,
    c1 = 133.5077,
    c2 = 137239,
    w_ref = 344010,
    z = 1,
    a_0 = 3.72,
    MM = 0.075718) "NiOH+";

  //406
  constant DataRecord NiO_aq(
    R = Modelica.Constants.R/NiO_aq.MM,
    G_ref = -164599,
    H_ref = -235266,
    S_ref = 104.6,
    a1 = -5.9973e-06,
    a2 = -4718.97,
    a3 = 0.00042573,
    a4 = -96759,
    c1 = 49.1649,
    c2 = -41735,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.07471) "NiO,aq";

  //407
  constant DataRecord HNiO2_n1(
    R = Modelica.Constants.R/HNiO2_n1.MM,
    G_ref = -343004,
    H_ref = -496892,
    S_ref = -150.624,
    a1 = -3.407e-06,
    a2 = -4084.38,
    a3 = 0.00040036,
    a4 = -99383,
    c1 = 227.0088,
    c2 = 282127,
    w_ref = 906420,
    z = -1,
    a_0 = 3.72,
    MM = 0.091718) "HNiO2-";

  //408
  constant DataRecord NiO2_n2(
    R = Modelica.Constants.R/NiO2_n2.MM,
    G_ref = -268613,
    H_ref = -425931,
    S_ref = -162.758,
    a1 = -2.0142e-06,
    a2 = -3744.39,
    a3 = 0.000387,
    a4 = -100788,
    c1 = 129.3333,
    c2 = -276115,
    w_ref = 1589590,
    z = -2,
    a_0 = 3.72,
    MM = 0.09071) "NiO2-2";

  //409
  constant DataRecord o_TOLUATE_aq(
    R = Modelica.Constants.R/o_TOLUATE_aq.MM,
    G_ref = -206857,
    H_ref = -393589,
    S_ref = 146.858,
    a1 = 6.6259e-05,
    a2 = 12204.52,
    a3 = -8.4985e-05,
    a4 = -166724,
    c1 = 433.8975,
    c2 = -127470,
    w_ref = 459280,
    z = -1,
    a_0 = 3.72,
    MM = 0.13514) "o-TOLUATE,aq";

  //410
  constant DataRecord o_TOLUIC_ACID_aq(
    R = Modelica.Constants.R/o_TOLUIC_ACID_aq.MM,
    G_ref = -229158,
    H_ref = -387606,
    S_ref = 242.254,
    a1 = 7.1107e-05,
    a2 = 13302.9,
    a3 = -0.00010988,
    a4 = -171264,
    c1 = 430.2612,
    c2 = 281872,
    w_ref = -85140,
    z = 0,
    a_0 = 3.72,
    MM = 0.13614) "o-TOLUIC-ACID,aq";

  //411
  constant DataRecord O2_aq(
    R = Modelica.Constants.R/O2_aq.MM,
    G_ref = 16544,
    H_ref = -12134,
    S_ref = 108.951,
    a1 = 2.4221e-05,
    a2 = 2658.35,
    a3 = 0.0001361,
    a4 = -127265,
    c1 = 147.917,
    c2 = 350310,
    w_ref = -164980,
    z = 0,
    a_0 = 3.72,
    MM = 0.032) "O2,aq";

  //412
  constant DataRecord OCTANE_aq(
    R = Modelica.Constants.R/OCTANE_aq.MM,
    G_ref = 35899,
    H_ref = -248571,
    S_ref = 266.939,
    a1 = 8.9609e-05,
    a2 = 18624.07,
    a3 = -0.00049136,
    a4 = -193263,
    c1 = 522.1314,
    c2 = 1727595,
    w_ref = -404260,
    z = 0,
    a_0 = 3.72,
    MM = 0.11422) "OCTANE,aq";

  //413
  constant DataRecord OCTANOATE_aq(
    R = Modelica.Constants.R/OCTANOATE_aq.MM,
    G_ref = -321206,
    H_ref = -634211,
    S_ref = 245.601,
    a1 = 8.5049e-05,
    a2 = 16462.74,
    a3 = -0.00018167,
    a4 = -184326,
    c1 = 585.8027,
    c2 = -106880,
    w_ref = 309700,
    z = -1,
    a_0 = 3.72,
    MM = 0.1432) "OCTANOATE,aq";

  //414
  constant DataRecord OCTANOIC_ACID_aq(
    R = Modelica.Constants.R/OCTANOIC_ACID_aq.MM,
    G_ref = -349155,
    H_ref = -631993,
    S_ref = 346.854,
    a1 = 9.1723e-05,
    a2 = 17975.34,
    a3 = -0.00021595,
    a4 = -190581,
    c1 = 600.5048,
    c2 = 484428,
    w_ref = -15820,
    z = 0,
    a_0 = 3.72,
    MM = 0.14421) "OCTANOIC-ACID,aq";

  //415
  constant DataRecord OH_n1(
    R = Modelica.Constants.R/OH_n1.MM,
    G_ref = -157297,
    H_ref = -230024,
    S_ref = -10.711,
    a1 = 5.2413e-06,
    a2 = 30.88,
    a3 = 7.7082e-05,
    a4 = -116403,
    c1 = 17.3636,
    c2 = -432877,
    w_ref = 721570,
    z = -1,
    a_0 = 3.5,
    MM = 0.017008) "OH-";

  //416
  constant DataRecord OXALATE_aq(
    R = Modelica.Constants.R/OXALATE_aq.MM,
    G_ref = -674042,
    H_ref = -825085,
    S_ref = 45.606,
    a1 = 2.9043e-05,
    a2 = 4046.6,
    a3 = 3.6271e-05,
    a4 = -132997,
    c1 = -45.0274,
    c2 = -781517,
    w_ref = 1275280,
    z = -2,
    a_0 = 3.72,
    MM = 0.08802) "OXALATE,aq";

  //417
  constant DataRecord OXALIC_ACID_aq(
    R = Modelica.Constants.R/OXALIC_ACID_aq.MM,
    G_ref = -705590,
    H_ref = -814123,
    S_ref = 184.096,
    a1 = 3.5267e-05,
    a2 = 5179.5,
    a3 = 7.4751e-05,
    a4 = -137683,
    c1 = 106.6895,
    c2 = -113725,
    w_ref = -123720,
    z = 0,
    a_0 = 3.72,
    MM = 0.090036) "OXALIC-ACID,aq";

  //418
  constant DataRecord p_TOLUATE_aq(
    R = Modelica.Constants.R/p_TOLUATE_aq.MM,
    G_ref = -219576,
    H_ref = -402333,
    S_ref = 160.247,
    a1 = 6.6191e-05,
    a2 = 12188.45,
    a3 = -8.4471e-05,
    a4 = -166657,
    c1 = 342.6683,
    c2 = -138574,
    w_ref = 439150,
    z = -1,
    a_0 = 3.72,
    MM = 0.13514) "p-TOLUATE,aq";

  //419
  constant DataRecord p_TOLUIC_ACID_aq(
    R = Modelica.Constants.R/p_TOLUIC_ACID_aq.MM,
    G_ref = -244513,
    H_ref = -402459,
    S_ref = 243.927,
    a1 = 7.111e-05,
    a2 = 13301.94,
    a3 = -0.00010944,
    a4 = -171259,
    c1 = 381.5419,
    c2 = 221518,
    w_ref = -84060,
    z = 0,
    a_0 = 3.72,
    MM = 0.13614) "p-TOLUIC-ACID,aq";

  //420
  constant DataRecord PENTANE_aq(
    R = Modelica.Constants.R/PENTANE_aq.MM,
    G_ref = 8912,
    H_ref = -173887,
    S_ref = 198.74,
    a1 = 6.282e-05,
    a2 = 12082.18,
    a3 = -0.00023409,
    a4 = -166218,
    c1 = 373.2421,
    c2 = 1177022,
    w_ref = -300960,
    z = 0,
    a_0 = 3.72,
    MM = 0.072145) "PENTANE,aq";

  //421
  constant DataRecord PENTANOATE_aq(
    R = Modelica.Constants.R/PENTANOATE_aq.MM,
    G_ref = -345724,
    H_ref = -562246,
    S_ref = 160.247,
    a1 = 5.8285e-05,
    a2 = 10395.4,
    a3 = -4.3526e-05,
    a4 = -159243,
    c1 = 362.2574,
    c2 = -136143,
    w_ref = 439150,
    z = -1,
    a_0 = 3.72,
    MM = 0.10112) "PENTANOATE,aq";

  //422
  constant DataRecord PENTANOIC_ACID_aq(
    R = Modelica.Constants.R/PENTANOIC_ACID_aq.MM,
    G_ref = -373380,
    H_ref = -559359,
    S_ref = 262.755,
    a1 = 6.4683e-05,
    a2 = 11847.25,
    a3 = -7.6877e-05,
    a4 = -165247,
    c1 = 381.0197,
    c2 = 219451,
    w_ref = -71550,
    z = 0,
    a_0 = 3.72,
    MM = 0.10213) "PENTANOIC-ACID,aq";

  //423
  constant DataRecord PHENOL_aq(
    R = Modelica.Constants.R/PHENOL_aq.MM,
    G_ref = -51108,
    H_ref = -151963,
    S_ref = 190.372,
    a1 = 5.6316e-05,
    a2 = 9950.22,
    a3 = -3.3271e-05,
    a4 = -157402,
    c1 = 283.4723,
    c2 = 99161,
    w_ref = -119540,
    z = 0,
    a_0 = 3.72,
    MM = 0.094107) "PHENOL,aq";

  //424
  constant DataRecord PHENYLALANINE_aq(
    R = Modelica.Constants.R/PHENYLALANINE_aq.MM,
    G_ref = -207108,
    H_ref = -460575,
    S_ref = 221.334,
    a1 = 7.6537e-05,
    a2 = 11900.97,
    a3 = 0.00053033,
    a4 = -165469,
    c1 = 348.0226,
    c2 = 90027,
    w_ref = -192920,
    z = 0,
    a_0 = 3.72,
    MM = 0.16518) "PHENYLALANINE,aq";

  //425
  constant DataRecord PIMELATE_aq(
    R = Modelica.Constants.R/PIMELATE_aq.MM,
    G_ref = -664670,
    H_ref = -983073,
    S_ref = 165.686,
    a1 = 7.066e-05,
    a2 = 13200.56,
    a3 = -0.00010728,
    a4 = -170841,
    c1 = 112.027,
    c2 = -174732,
    w_ref = 1093990,
    z = -2,
    a_0 = 3.72,
    MM = 0.15815) "PIMELATE,aq";

  //426
  constant DataRecord PIMELIC_ACID_aq(
    R = Modelica.Constants.R/PIMELIC_ACID_aq.MM,
    G_ref = -721238,
    H_ref = -978173,
    S_ref = 372.794,
    a1 = 8.2924e-05,
    a2 = 15980.03,
    a3 = -0.0001704,
    a4 = -182330,
    c1 = 370.0166,
    c2 = 197543,
    w_ref = 1380,
    z = 0,
    a_0 = 3.72,
    MM = 0.16016) "PIMELIC-ACID,aq";

  //427
  constant DataRecord PO4_n3(
    R = Modelica.Constants.R/PO4_n3.MM,
    G_ref = -1018804,
    H_ref = -1277375,
    S_ref = -221.752,
    a1 = -2.1999e-06,
    a2 = -3789.7,
    a3 = 0.00038881,
    a4 = -100604,
    c1 = -63.429,
    c2 = -1188905,
    w_ref = 2347810,
    z = -3,
    a_0 = 4,
    MM = 0.09497) "PO4-3";

  //428
  constant DataRecord PROPANE_aq(
    R = Modelica.Constants.R/PROPANE_aq.MM,
    G_ref = -8213,
    H_ref = -127570,
    S_ref = 139.536,
    a1 = 4.5053e-05,
    a2 = 7396.68,
    a3 = -2.4594e-05,
    a4 = -146846,
    c1 = 353.7405,
    c2 = 193853,
    w_ref = -69040,
    z = 0,
    a_0 = 3.72,
    MM = 0.044093) "PROPANE,aq";

  //429
  constant DataRecord PROPANAMINE_aq(
    R = Modelica.Constants.R/PROPANAMINE_aq.MM,
    G_ref = 29330,
    H_ref = -128365,
    S_ref = 170.707,
    a1 = 4.9446e-05,
    a2 = 6485.74,
    a3 = 0.00048582,
    a4 = -143080,
    c1 = 305.2784,
    c2 = 50995,
    w_ref = -116270,
    z = 0,
    a_0 = 3.72,
    MM = 0.059108) "PROPANAMINE,aq";

  //430
  constant DataRecord PROPANOATE_aq(
    R = Modelica.Constants.R/PROPANOATE_aq.MM,
    G_ref = -363088,
    H_ref = -513084,
    S_ref = 110.876,
    a1 = 4.0581e-05,
    a2 = 5077.03,
    a3 = 0.00037912,
    a4 = -137256,
    c1 = 218.8232,
    c2 = -175728,
    w_ref = 513630,
    z = -1,
    a_0 = 3.72,
    MM = 0.07307) "PROPANOATE,aq";

  //431
  constant DataRecord PROPANOIC_ACID_aq(
    R = Modelica.Constants.R/PROPANOIC_ACID_aq.MM,
    G_ref = -390995,
    H_ref = -512414,
    S_ref = 206.69,
    a1 = 4.6048e-05,
    a2 = 7827.3,
    a3 = -3.2602e-05,
    a4 = -148628,
    c1 = 263.4665,
    c2 = -49790,
    w_ref = -62760,
    z = 0,
    a_0 = 3.72,
    MM = 0.074077) "PROPANOIC-ACID,aq";

  //432
  constant DataRecord PROPANOL_aq(
    R = Modelica.Constants.R/PROPANOL_aq.MM,
    G_ref = -175351,
    H_ref = -315139,
    S_ref = 173.218,
    a1 = 4.7457e-05,
    a2 = 5993.66,
    a3 = 0.0005066,
    a4 = -141047,
    c1 = 327.6666,
    c2 = 68814,
    w_ref = -120040,
    z = 0,
    a_0 = 3.72,
    MM = 0.060093) "PROPANOL,aq";

  //433
  constant DataRecord Pb_p2(
    R = Modelica.Constants.R/Pb_p2.MM,
    G_ref = -23891,
    H_ref = 920,
    S_ref = 17.573,
    a1 = -2.13e-08,
    a2 = -3260.97,
    a3 = 0.00036875,
    a4 = -102793,
    c1 = 36.2435,
    c2 = -235208,
    w_ref = 451370,
    z = 2,
    a_0 = 4.5,
    MM = 0.20719) "Pb+2";

  //434
  constant DataRecord PbCl_p1(
    R = Modelica.Constants.R/PbCl_p1.MM,
    G_ref = -163385,
    H_ref = -161628,
    S_ref = 117.152,
    a1 = 1.245e-05,
    a2 = -214.64,
    a3 = 0.00024872,
    a4 = -115382,
    c1 = 42.7529,
    c2 = -85203,
    w_ref = 53600,
    z = 1,
    a_0 = 3.72,
    MM = 0.24264) "PbCl+";

  //435
  constant DataRecord PbCl2_aq(
    R = Modelica.Constants.R/PbCl2_aq.MM,
    G_ref = -297901,
    H_ref = -325097,
    S_ref = 196.648,
    a1 = 2.7794e-05,
    a2 = 3531.97,
    a3 = 0.00010147,
    a4 = -130871,
    c1 = 29.2403,
    c2 = -109922,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.27809) "PbCl2,aq";

  //436
  constant DataRecord PbCl3_n1(
    R = Modelica.Constants.R/PbCl3_n1.MM,
    G_ref = -427396,
    H_ref = -492457,
    S_ref = 246.856,
    a1 = 4.6253e-05,
    a2 = 8039.18,
    a3 = -7.5684e-05,
    a4 = -149507,
    c1 = 24.9776,
    c2 = -228388,
    w_ref = 307780,
    z = -1,
    a_0 = 3.72,
    MM = 0.31354) "PbCl3-";

  //437
  constant DataRecord PbCl4_n2(
    R = Modelica.Constants.R/PbCl4_n2.MM,
    G_ref = -557213,
    H_ref = -666130,
    S_ref = 276.144,
    a1 = 6.7687e-05,
    a2 = 13272.9,
    a3 = -0.0002814,
    a4 = -171142,
    c1 = 24.8258,
    c2 = -426797,
    w_ref = 925750,
    z = -2,
    a_0 = 3.72,
    MM = 0.34899) "PbCl4-2";

  //438
  constant DataRecord PbOH_p1(
    R = Modelica.Constants.R/PbOH_p1.MM,
    G_ref = -225727,
    H_ref = -288278,
    S_ref = -41.84,
    a1 = 1.0321e-05,
    a2 = -736.38,
    a3 = 0.00026971,
    a4 = -113223,
    c1 = 3.753,
    c2 = -297420,
    w_ref = 293010,
    z = 1,
    a_0 = 3.72,
    MM = 0.2242) "PbOH+";

  //439
  constant DataRecord PbO_aq(
    R = Modelica.Constants.R/PbO_aq.MM,
    G_ref = -164640,
    H_ref = -187025,
    S_ref = 92.048,
    a1 = 1.2094e-05,
    a2 = -302.75,
    a3 = 0.00025251,
    a4 = -115018,
    c1 = -78.3429,
    c2 = -484921,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22319) "PbO,aq";

  //440
  constant DataRecord HPbO2_n1(
    R = Modelica.Constants.R/HPbO2_n1.MM,
    G_ref = -338904,
    H_ref = -435973,
    S_ref = 75.312,
    a1 = 1.8578e-05,
    a2 = 1280.26,
    a3 = 0.0001903,
    a4 = -121562,
    c1 = -93.5542,
    c2 = -723560,
    w_ref = 567600,
    z = -1,
    a_0 = 3.72,
    MM = 0.2402) "HPbO2-";

  //441
  constant DataRecord PdOH_p1(
    R = Modelica.Constants.R/PdOH_p1.MM,
    G_ref = -54936,
    H_ref = -100918,
    S_ref = -14.226,
    a1 = 5.13e-07,
    a2 = -3127.16,
    a3 = 0.00036274,
    a4 = -103341,
    c1 = 80.5575,
    c2 = -17870,
    w_ref = 253680,
    z = 1,
    a_0 = 3.72,
    MM = 0.12341) "PdOH+";

  //442
  constant DataRecord PdSO4_aq(
    R = Modelica.Constants.R/PdSO4_aq.MM,
    G_ref = -581576,
    H_ref = -725715,
    S_ref = -3.766,
    a1 = 9.6324e-06,
    a2 = -903.12,
    a3 = 0.00027594,
    a4 = -112537,
    c1 = -24.5442,
    c2 = -297930,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.20246) "PdSO4,aq";

  //443
  constant DataRecord Pd_SO42_n2(
    R = Modelica.Constants.R/Pd_SO42_n2.MM,
    G_ref = -1337374,
    H_ref = -1627158,
    S_ref = 80.333,
    a1 = 2.7204e-05,
    a2 = 3385.53,
    a3 = 0.00010778,
    a4 = -130265,
    c1 = 54.087,
    c2 = -420149,
    w_ref = 1222610,
    z = -2,
    a_0 = 3.72,
    MM = 0.29852) "Pd(SO4)2-2";

  //444
  constant DataRecord Pd_SO43_n4(
    R = Modelica.Constants.R/Pd_SO43_n4.MM,
    G_ref = -2091582,
    H_ref = -2526927,
    S_ref = 164.85,
    a1 = 4.48e-05,
    a2 = 7683.75,
    a3 = -6.1555e-05,
    a4 = -148034,
    c1 = 134.6867,
    c2 = -537761,
    w_ref = 2464790,
    z = -4,
    a_0 = 3.72,
    MM = 0.39458) "Pd(SO4)3-4";

  //445
  constant DataRecord Pd_p2(
    R = Modelica.Constants.R/Pd_p2.MM,
    G_ref = 176565,
    H_ref = 177987,
    S_ref = -88.282,
    a1 = -1.5192e-06,
    a2 = -3623.43,
    a3 = 0.00038225,
    a4 = -101290,
    c1 = 71.9259,
    c2 = -175372,
    w_ref = 651830,
    z = 2,
    a_0 = 3.72,
    MM = 0) "Pd+2";

  //446
  constant DataRecord PdCl_p1(
    R = Modelica.Constants.R/PdCl_p1.MM,
    G_ref = 10460,
    H_ref = -20585,
    S_ref = -16.318,
    a1 = 1.1183e-05,
    a2 = -526.22,
    a3 = 0.00026151,
    a4 = -114093,
    c1 = 101.057,
    c2 = 52304,
    w_ref = 257020,
    z = 1,
    a_0 = 3.72,
    MM = 0.14185) "PdCl+";

  //447
  constant DataRecord PdCl2_aq(
    R = Modelica.Constants.R/PdCl2_aq.MM,
    G_ref = -147109,
    H_ref = -219493,
    S_ref = 17.573,
    a1 = 2.5913e-05,
    a2 = 3071.22,
    a3 = 0.00011993,
    a4 = -128968,
    c1 = 117.3315,
    c2 = 195184,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1773) "PdCl2,aq";

  //448
  constant DataRecord PdCl3_n1(
    R = Modelica.Constants.R/PdCl3_n1.MM,
    G_ref = -292085,
    H_ref = -419237,
    S_ref = 10.46,
    a1 = 4.5629e-05,
    a2 = 7886.71,
    a3 = -6.9643e-05,
    a4 = -148871,
    c1 = 196.4187,
    c2 = 253094,
    w_ref = 665050,
    z = -1,
    a_0 = 3.72,
    MM = 0.21275) "PdCl3-";

  //449
  constant DataRecord PdCl4_n2(
    R = Modelica.Constants.R/PdCl4_n2.MM,
    G_ref = -434718,
    H_ref = -632453,
    S_ref = -48.953,
    a1 = 6.7593e-05,
    a2 = 13249.14,
    a3 = -0.00028024,
    a4 = -171042,
    c1 = 258.0879,
    c2 = 226036,
    w_ref = 1418960,
    z = -2,
    a_0 = 3.72,
    MM = 0.2482) "PdCl4-2";

  //450
  constant DataRecord PdO_aq(
    R = Modelica.Constants.R/PdO_aq.MM,
    G_ref = -48116,
    H_ref = -101755,
    S_ref = -39.748,
    a1 = -1.9895e-06,
    a2 = -3740.41,
    a3 = 0.00038727,
    a4 = -100805,
    c1 = 13.1194,
    c2 = -167021,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1224) "PdO,aq";

  //451
  constant DataRecord Pr_p2(
    R = Modelica.Constants.R/Pr_p2.MM,
    G_ref = -387857,
    H_ref = -369866,
    S_ref = 2.929,
    a1 = 4.452e-07,
    a2 = -3143.36,
    a3 = 0.00036332,
    a4 = -103274,
    c1 = 36.6313,
    c2 = -229237,
    w_ref = 436980,
    z = 2,
    a_0 = 3.72,
    MM = 0.14091) "Pr+2";

  //452
  constant DataRecord Pr_p3(
    R = Modelica.Constants.R/Pr_p3.MM,
    G_ref = -680318,
    H_ref = -706259,
    S_ref = -209.2,
    a1 = -1.3414e-05,
    a2 = -5978.68,
    a3 = 0.00035701,
    a4 = -91554,
    c1 = 16.5716,
    c2 = -472139,
    w_ref = 977760,
    z = 3,
    a_0 = 9,
    MM = 0.14091) "Pr+3";

  //453
  constant DataRecord Pr_p4(
    R = Modelica.Constants.R/Pr_p4.MM,
    G_ref = -304177,
    H_ref = -371121,
    S_ref = -412.124,
    a1 = -1.7883e-05,
    a2 = -7620.24,
    a3 = 0.00053956,
    a4 = -84768,
    c1 = 168.7646,
    c2 = -121851,
    w_ref = 1535820,
    z = 4,
    a_0 = 3.72,
    MM = 0.14091) "Pr+4";

  //454
  constant DataRecord PtOH_p1(
    R = Modelica.Constants.R/PtOH_p1.MM,
    G_ref = 6527,
    H_ref = -44769,
    S_ref = -28.033,
    a1 = -2.757e-07,
    a2 = -3320,
    a3 = 0.00037036,
    a4 = -102546,
    c1 = 91.7823,
    c2 = 14510,
    w_ref = 274340,
    z = 1,
    a_0 = 3.72,
    MM = 0.2121) "PtOH+";

  //455
  constant DataRecord PtSO4_aq(
    R = Modelica.Constants.R/PtSO4_aq.MM,
    G_ref = -503335,
    H_ref = -652076,
    S_ref = -15.062,
    a1 = 8.4739e-06,
    a2 = -1187,
    a3 = 0.00028727,
    a4 = -111361,
    c1 = -27.8056,
    c2 = -309265,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.29115) "PtSO4,aq";

  //456
  constant DataRecord Pt_SO42_n2(
    R = Modelica.Constants.R/Pt_SO42_n2.MM,
    G_ref = -1258003,
    H_ref = -1552599,
    S_ref = 68.618,
    a1 = 2.6263e-05,
    a2 = 3157,
    a3 = 0.00011645,
    a4 = -129319,
    c1 = 52.3435,
    c2 = -431822,
    w_ref = 1240180,
    z = -2,
    a_0 = 3.72,
    MM = 0.38721) "Pt(SO4)2-2";

  //457
  constant DataRecord Pt_SO43_n4(
    R = Modelica.Constants.R/Pt_SO43_n4.MM,
    G_ref = -2009910,
    H_ref = -2462828,
    S_ref = 151.879,
    a1 = 4.4017e-05,
    a2 = 7491.87,
    a3 = -5.3873e-05,
    a4 = -147239,
    c1 = 131.5445,
    c2 = -554384,
    w_ref = 2482580,
    z = -4,
    a_0 = 3.72,
    MM = 0.48327) "Pt(SO4)3-4";

  //458
  constant DataRecord Pt_p2(
    R = Modelica.Constants.R/Pt_p2.MM,
    G_ref = 257734,
    H_ref = 254806,
    S_ref = -98.742,
    a1 = -2.7815e-06,
    a2 = -3932.33,
    a3 = 0.00039452,
    a4 = -100014,
    c1 = 70.3029,
    c2 = -186707,
    w_ref = 669650,
    z = 2,
    a_0 = 3.72,
    MM = 0.19509) "Pt+2";

  //459
  constant DataRecord PtCl_p1(
    R = Modelica.Constants.R/PtCl_p1.MM,
    G_ref = 76902,
    H_ref = 41966,
    S_ref = -29.288,
    a1 = 9.7751e-06,
    a2 = -869.98,
    a3 = 0.00027496,
    a4 = -112671,
    c1 = 95.4831,
    c2 = 27376,
    w_ref = 274340,
    z = 1,
    a_0 = 3.72,
    MM = 0.23054) "PtCl+";

  //460
  constant DataRecord PtCl2_aq(
    R = Modelica.Constants.R/PtCl2_aq.MM,
    G_ref = -93261,
    H_ref = -172004,
    S_ref = 0.418,
    a1 = 2.4278e-05,
    a2 = 2671.86,
    a3 = 0.00013566,
    a4 = -127315,
    c1 = 103.3373,
    c2 = 146545,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.26599) "PtCl2,aq";

  //461
  constant DataRecord PtCl3_n1(
    R = Modelica.Constants.R/PtCl3_n1.MM,
    G_ref = -241668,
    H_ref = -377188,
    S_ref = -12.97,
    a1 = 4.393e-05,
    a2 = 7470.41,
    a3 = -5.2961e-05,
    a4 = -147151,
    c1 = 176.0949,
    c2 = 170883,
    w_ref = 701200,
    z = -1,
    a_0 = 3.72,
    MM = 0.30144) "PtCl3-";

  //462
  constant DataRecord PtCl4_n2(
    R = Modelica.Constants.R/PtCl4_n2.MM,
    G_ref = -381623,
    H_ref = -590613,
    S_ref = -82.843,
    a1 = 6.5732e-05,
    a2 = 12794.34,
    a3 = -0.00026234,
    a4 = -169159,
    c1 = 226.5142,
    c2 = 100391,
    w_ref = 1468580,
    z = -2,
    a_0 = 3.72,
    MM = 0.33689) "PtCl4-2";

  //463
  constant DataRecord PtO_aq(
    R = Modelica.Constants.R/PtO_aq.MM,
    G_ref = -4728,
    H_ref = -64015,
    S_ref = -54.392,
    a1 = -2.9058e-06,
    a2 = -3962.46,
    a3 = 0.00039567,
    a4 = -99889,
    c1 = 22.9279,
    c2 = -132930,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.21109) "PtO,aq";

  //464
  constant DataRecord Ra_p2(
    R = Modelica.Constants.R/Ra_p2.MM,
    G_ref = -561493,
    H_ref = -527602,
    S_ref = 53.974,
    a1 = 1.5648e-06,
    a2 = -2870.14,
    a3 = 0.00035262,
    a4 = -104403,
    c1 = 24.0563,
    c2 = -248839,
    w_ref = 361660,
    z = 2,
    a_0 = 5,
    MM = 0.226) "Ra+2";

  //465
  constant DataRecord Rb_p1(
    R = Modelica.Constants.R/Rb_p1.MM,
    G_ref = -283675,
    H_ref = -251124,
    S_ref = 120.499,
    a1 = 1.7955e-05,
    a2 = -378.28,
    a3 = 0.00030991,
    a4 = -114709,
    c1 = 24.235,
    c2 = -152536,
    w_ref = 62840,
    z = 1,
    a_0 = 2.5,
    MM = 0.08547) "Rb+";

  //466
  constant DataRecord RbBr_aq(
    R = Modelica.Constants.R/RbBr_aq.MM,
    G_ref = -380786,
    H_ref = -358694,
    S_ref = 226.773,
    a1 = 3.5425e-05,
    a2 = 5395.35,
    a3 = 2.8229e-05,
    a4 = -138574,
    c1 = -21.6396,
    c2 = -290520,
    w_ref = -4180,
    z = 0,
    a_0 = 3.72,
    MM = 0.16538) "RbBr,aq";

  //467
  constant DataRecord RbCl_aq(
    R = Modelica.Constants.R/RbCl_aq.MM,
    G_ref = -409488,
    H_ref = -405011,
    S_ref = 203.175,
    a1 = 3.1188e-05,
    a2 = 4360.86,
    a3 = 6.889e-05,
    a4 = -134298,
    c1 = -23.6262,
    c2 = -297424,
    w_ref = -4180,
    z = 0,
    a_0 = 3.72,
    MM = 0.12092) "RbCl,aq";

  //468
  constant DataRecord RbF_aq(
    R = Modelica.Constants.R/RbF_aq.MM,
    G_ref = -570907,
    H_ref = -584547,
    S_ref = 132.214,
    a1 = 1.8262e-05,
    a2 = 1204.49,
    a3 = 0.00019295,
    a4 = -121252,
    c1 = -13.4708,
    c2 = -263333,
    w_ref = -420,
    z = 0,
    a_0 = 3.72,
    MM = 0.10447) "RbF,aq";

  //469
  constant DataRecord RbI_aq(
    R = Modelica.Constants.R/RbI_aq.MM,
    G_ref = -330118,
    H_ref = -300076,
    S_ref = 235.559,
    a1 = 4.4026e-05,
    a2 = 7495.43,
    a3 = -5.4317e-05,
    a4 = -147256,
    c1 = -15.972,
    c2 = -272023,
    w_ref = -420,
    z = 0,
    a_0 = 3.72,
    MM = 0.21237) "RbI,aq";

  //470
  constant DataRecord RbOH_aq(
    R = Modelica.Constants.R/RbOH_aq.MM,
    G_ref = -439738,
    H_ref = -472792,
    S_ref = 133.888,
    a1 = 1.9194e-05,
    a2 = 1429.59,
    a3 = 0.00018465,
    a4 = -122181,
    c1 = -42.7873,
    c2 = -361343,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.10248) "RbOH,aq";

  //471
  constant DataRecord ReO4_n1(
    R = Modelica.Constants.R/ReO4_n1.MM,
    G_ref = -694544,
    H_ref = -787429,
    S_ref = 201.25,
    a1 = 3.6026e-05,
    a2 = 8392.31,
    a3 = -0.00070099,
    a4 = -150963,
    c1 = 79.2031,
    c2 = -296227,
    w_ref = 376730,
    z = -1,
    a_0 = 3.72,
    MM = 0.2502) "ReO4-";

  //472
  constant DataRecord RhOH_p1(
    R = Modelica.Constants.R/RhOH_p1.MM,
    G_ref = -76525,
    H_ref = -145478,
    S_ref = -52.718,
    a1 = -1.7046e-06,
    a2 = -3668.87,
    a3 = 0.00038405,
    a4 = -101102,
    c1 = 111.6433,
    c2 = 72467,
    w_ref = 308950,
    z = 1,
    a_0 = 3.72,
    MM = 0.11992) "RhOH+";

  //473
  constant DataRecord RhOH_p2(
    R = Modelica.Constants.R/RhOH_p2.MM,
    G_ref = -3431,
    H_ref = -61170,
    S_ref = -135.98,
    a1 = 9.4014e-06,
    a2 = -960.06,
    a3 = 0.00027825,
    a4 = -112299,
    c1 = 85.2072,
    c2 = -127817,
    w_ref = 647470,
    z = 2,
    a_0 = 3.72,
    MM = 0.11992) "RhOH+2";

  //474
  constant DataRecord RhSO4_aq(
    R = Modelica.Constants.R/RhSO4_aq.MM,
    G_ref = -642244,
    H_ref = -794123,
    S_ref = -35.564,
    a1 = 6.3576e-06,
    a2 = -1700.38,
    a3 = 0.00030675,
    a4 = -109240,
    c1 = -33.6419,
    c2 = -329549,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19897) "RhSO4,aq";

  //475
  constant DataRecord RhSO4_p1(
    R = Modelica.Constants.R/RhSO4_p1.MM,
    G_ref = -533460,
    H_ref = -698728,
    S_ref = -146.858,
    a1 = -5.289e-07,
    a2 = -3381.63,
    a3 = 0.00037276,
    a4 = -102290,
    c1 = -41.5358,
    c2 = -506908,
    w_ref = 455640,
    z = 1,
    a_0 = 3.72,
    MM = 0.19897) "RhSO4+";

  //476
  constant DataRecord Rh_SO42_n1(
    R = Modelica.Constants.R/Rh_SO42_n1.MM,
    G_ref = -1280722,
    H_ref = -1571510,
    S_ref = 5.858,
    a1 = 1.5659e-05,
    a2 = 567.35,
    a3 = 0.00021836,
    a4 = -118616,
    c1 = -56.72,
    c2 = -629470,
    w_ref = 673540,
    z = -1,
    a_0 = 3.72,
    MM = 0.29503) "Rh(SO4)2-";

  //477
  constant DataRecord Rh_SO42_n2(
    R = Modelica.Constants.R/Rh_SO42_n2.MM,
    G_ref = -1397874,
    H_ref = -1695775,
    S_ref = 46.861,
    a1 = 2.7409e-05,
    a2 = 3436.49,
    a3 = 0.00010553,
    a4 = -130474,
    c1 = 49.5206,
    c2 = -452110,
    w_ref = 1272860,
    z = -2,
    a_0 = 3.72,
    MM = 0.29503) "Rh(SO4)2-2";

  //478
  constant DataRecord Rh_SO43_n3(
    R = Modelica.Constants.R/Rh_SO43_n3.MM,
    G_ref = -2023801,
    H_ref = -2439690,
    S_ref = 158.574,
    a1 = 3.4826e-05,
    a2 = 5247.41,
    a3 = 3.4443e-05,
    a4 = -137963,
    c1 = 9.3086,
    c2 = -752028,
    w_ref = 1773010,
    z = -3,
    a_0 = 3.72,
    MM = 0.39109) "Rh(SO4)3-3";

  //479
  constant DataRecord Rh_SO43_n4(
    R = Modelica.Constants.R/Rh_SO43_n4.MM,
    G_ref = -2150994,
    H_ref = -2595335,
    S_ref = 128.867,
    a1 = 4.2603e-05,
    a2 = 7147.15,
    a3 = -4.1706e-05,
    a4 = -145817,
    c1 = 129.1011,
    c2 = -574668,
    w_ref = 2519400,
    z = -4,
    a_0 = 3.72,
    MM = 0.39109) "Rh(SO4)3-4";

  //480
  constant DataRecord Rh_p2(
    R = Modelica.Constants.R/Rh_p2.MM,
    G_ref = 115897,
    H_ref = 110458,
    S_ref = -117.57,
    a1 = -5.3551e-06,
    a2 = -4560.89,
    a3 = 0.00041925,
    a4 = -97416,
    c1 = 59.7028,
    c2 = -206991,
    w_ref = 617930,
    z = 2,
    a_0 = 3.72,
    MM = 0) "Rh+2";

  //481
  constant DataRecord Rh_p3(
    R = Modelica.Constants.R/Rh_p3.MM,
    G_ref = 219451,
    H_ref = 179201,
    S_ref = -299.574,
    a1 = -1.3751e-05,
    a2 = -6611.56,
    a3 = 0.00049996,
    a4 = -88935,
    c1 = 54.4887,
    c2 = -384351,
    w_ref = 1115200,
    z = 3,
    a_0 = 3.72,
    MM = 0) "Rh+3";

  //482
  constant DataRecord RhCl_p1(
    R = Modelica.Constants.R/RhCl_p1.MM,
    G_ref = -14142,
    H_ref = -53011,
    S_ref = -53.137,
    a1 = 7.2128e-06,
    a2 = -1493.65,
    a3 = 0.00029906,
    a4 = -110094,
    c1 = 85.4285,
    c2 = -18644,
    w_ref = 308950,
    z = 1,
    a_0 = 3.72,
    MM = 0.13836) "RhCl+";

  //483
  constant DataRecord RhCl_p2(
    R = Modelica.Constants.R/RhCl_p2.MM,
    G_ref = 76693,
    H_ref = 11724,
    S_ref = -205.434,
    a1 = -2.5167e-06,
    a2 = -3868.57,
    a3 = 0.00039221,
    a4 = -100278,
    c1 = 122.5753,
    c2 = -31430,
    w_ref = 752070,
    z = 2,
    a_0 = 3.72,
    MM = 0.13836) "RhCl+2";

  //484
  constant DataRecord RhCl2_aq(
    R = Modelica.Constants.R/RhCl2_aq.MM,
    G_ref = -142130,
    H_ref = -226982,
    S_ref = -30.543,
    a1 = 2.1293e-05,
    a2 = 1942.34,
    a3 = 0.00016447,
    a4 = -124298,
    c1 = 77.5023,
    c2 = 56748,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17381) "RhCl2,aq";

  //485
  constant DataRecord RhCl2_p1(
    R = Modelica.Constants.R/RhCl2_p1.MM,
    G_ref = -61840,
    H_ref = -166732,
    S_ref = -162.758,
    a1 = 1.0436e-05,
    a2 = -708.9,
    a3 = 0.0002687,
    a4 = -113340,
    c1 = 120.7912,
    c2 = 50681,
    w_ref = 476260,
    z = 1,
    a_0 = 3.72,
    MM = 0.17381) "RhCl2+";

  //486
  constant DataRecord RhCl3_aq(
    R = Modelica.Constants.R/RhCl3_aq.MM,
    G_ref = -193259,
    H_ref = -354343,
    S_ref = -174.473,
    a1 = 2.4254e-05,
    a2 = 2665.17,
    a3 = 0.00013609,
    a4 = -127286,
    c1 = 21.466,
    c2 = -138013,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.20926) "RhCl3,aq";

  //487
  constant DataRecord RhCl3_n1(
    R = Modelica.Constants.R/RhCl3_n1.MM,
    G_ref = -265810,
    H_ref = -410869,
    S_ref = -55.647,
    a1 = 4.0816e-05,
    a2 = 6709.38,
    a3 = -2.2916e-05,
    a4 = -144005,
    c1 = 138.2042,
    c2 = 19113,
    w_ref = 763870,
    z = -1,
    a_0 = 3.72,
    MM = 0.20926) "RhCl3-";

  //488
  constant DataRecord RhCl4_n1(
    R = Modelica.Constants.R/RhCl4_n1.MM,
    G_ref = -324260,
    H_ref = -561493,
    S_ref = -253.132,
    a1 = 4.5124e-05,
    a2 = 7763.37,
    a3 = -6.4806e-05,
    a4 = -148365,
    c1 = -11.711,
    c2 = -597521,
    w_ref = 1062320,
    z = -1,
    a_0 = 3.72,
    MM = 0.24471) "RhCl4-";

  //489
  constant DataRecord RhCl4_n2(
    R = Modelica.Constants.R/RhCl4_n2.MM,
    G_ref = -390158,
    H_ref = -614128,
    S_ref = -143.511,
    a1 = 6.2342e-05,
    a2 = 11967.03,
    a3 = -0.00022987,
    a4 = -165741,
    c1 = 168.4039,
    c2 = -131562,
    w_ref = 1562220,
    z = -2,
    a_0 = 3.72,
    MM = 0.24471) "RhCl4-2";

  //490
  constant DataRecord RhO_aq(
    R = Modelica.Constants.R/RhO_aq.MM,
    G_ref = -30208,
    H_ref = -94307,
    S_ref = -81.17,
    a1 = -4.5087e-06,
    a2 = -4354.46,
    a3 = 0.00041119,
    a4 = -98270,
    c1 = 40.8283,
    c2 = -70714,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.11891) "RhO,aq";

  //491
  constant DataRecord RhO_p1(
    R = Modelica.Constants.R/RhO_p1.MM,
    G_ref = 13138,
    H_ref = -41882,
    S_ref = -78.241,
    a1 = 9.6512e-06,
    a2 = -900.19,
    a3 = 0.00027615,
    a4 = -112550,
    c1 = -38.9317,
    c2 = -463612,
    w_ref = 348690,
    z = 1,
    a_0 = 3.72,
    MM = 0.11891) "RhO+";

  //492
  constant DataRecord Rn_aq(
    R = Modelica.Constants.R/Rn_aq.MM,
    G_ref = 11673,
    H_ref = -20920,
    S_ref = 66.944,
    a1 = 3.8017e-05,
    a2 = 6026.84,
    a3 = 3.699e-06,
    a4 = -141189,
    c1 = 219.9838,
    c2 = 580425,
    w_ref = -101380,
    z = 0,
    a_0 = 3.72,
    MM = 0.222) "Rn,aq";

  //493
  constant DataRecord RuOH_p1(
    R = Modelica.Constants.R/RuOH_p1.MM,
    G_ref = -43806,
    H_ref = -110081,
    S_ref = -44.35,
    a1 = -1.2301e-06,
    a2 = -3552.97,
    a3 = 0.00037949,
    a4 = -101579,
    c1 = 104.8933,
    c2 = 52865,
    w_ref = 296900,
    z = 1,
    a_0 = 3.72,
    MM = 0.11808) "RuOH+";

  //494
  constant DataRecord RuOH_p2(
    R = Modelica.Constants.R/RuOH_p2.MM,
    G_ref = -51003,
    H_ref = -107194,
    S_ref = -131.796,
    a1 = 9.501e-06,
    a2 = -936.55,
    a3 = 0.00027755,
    a4 = -112399,
    c1 = 81.8662,
    c2 = -138043,
    w_ref = 643160,
    z = 2,
    a_0 = 3.72,
    MM = 0.11808) "RuOH+2";

  //495
  constant DataRecord RuO4_n2(
    R = Modelica.Constants.R/RuO4_n2.MM,
    G_ref = -299574,
    H_ref = -461077,
    S_ref = 27.614,
    a1 = 2.908e-05,
    a2 = 3844.55,
    a3 = 8.9483e-05,
    a4 = -132164,
    c1 = 27.9102,
    c2 = -536912,
    w_ref = 1303110,
    z = -1,
    a_0 = 3.72,
    MM = 0.16507) "RuO4-2";

  //496
  constant DataRecord RuSO4_aq(
    R = Modelica.Constants.R/RuSO4_aq.MM,
    G_ref = -607517,
    H_ref = -756467,
    S_ref = -28.87,
    a1 = 7.063e-06,
    a2 = -1530.13,
    a3 = 0.00030048,
    a4 = -109943,
    c1 = -31.7043,
    c2 = -322817,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19713) "RuSO4,aq";

  //497
  constant DataRecord RuSO4_p1(
    R = Modelica.Constants.R/RuSO4_p1.MM,
    G_ref = -582413,
    H_ref = -746007,
    S_ref = -143.511,
    a1 = 1.5134e-06,
    a2 = -2882.78,
    a3 = 0.00035311,
    a4 = -104353,
    c1 = -38.4664,
    c2 = -494126,
    w_ref = 449030,
    z = 1,
    a_0 = 3.72,
    MM = 0.19713) "RuSO4+";

  //498
  constant DataRecord Ru_SO42_n1(
    R = Modelica.Constants.R/Ru_SO42_n1.MM,
    G_ref = -1330094,
    H_ref = -1618790,
    S_ref = 9.205,
    a1 = 1.7424e-05,
    a2 = 997.88,
    a3 = 0.00020151,
    a4 = -120395,
    c1 = -53.5682,
    c2 = -616684,
    w_ref = 667850,
    z = -1,
    a_0 = 3.72,
    MM = 0.29319) "Ru(SO4)2-";

  //499
  constant DataRecord Ru_SO42_n2(
    R = Modelica.Constants.R/Ru_SO42_n2.MM,
    G_ref = -1361474,
    H_ref = -1656446,
    S_ref = 53.974,
    a1 = 2.5121e-05,
    a2 = 2879.18,
    a3 = 0.00012721,
    a4 = -128173,
    c1 = 50.5712,
    c2 = -445374,
    w_ref = 1263230,
    z = -2,
    a_0 = 3.72,
    MM = 0.28334) "Ru(SO4)2-2";

  //500
  constant DataRecord Ru_SO43_n3(
    R = Modelica.Constants.R/Ru_SO43_n3.MM,
    G_ref = -2071917,
    H_ref = -2486133,
    S_ref = 161.921,
    a1 = 3.631e-05,
    a2 = 5610.7,
    a3 = 1.9887e-05,
    a4 = -139465,
    c1 = 12.4574,
    c2 = -739242,
    w_ref = 1767280,
    z = -3,
    a_0 = 3.72,
    MM = 0.38925) "Ru(SO4)3-3";

  //501
  constant DataRecord Ru_SO43_n4(
    R = Modelica.Constants.R/Ru_SO43_n4.MM,
    G_ref = -2112083,
    H_ref = -2553077,
    S_ref = 136.817,
    a1 = 4.3073e-05,
    a2 = 7261.62,
    a3 = -4.4878e-05,
    a4 = -146289,
    c1 = 129.9605,
    c2 = -567932,
    w_ref = 2507720,
    z = -4,
    a_0 = 3.72,
    MM = 0.38925) "Ru(SO4)3-4";

  //502
  constant DataRecord Ru_p2(
    R = Modelica.Constants.R/Ru_p2.MM,
    G_ref = 150206,
    H_ref = 147486,
    S_ref = -111.294,
    a1 = -4.5639e-06,
    a2 = -4368.72,
    a3 = 0.00041189,
    a4 = -98207,
    c1 = 60.8885,
    c2 = -200259,
    w_ref = 609780,
    z = 2,
    a_0 = 3.72,
    MM = 0) "Ru+2";

  //503
  constant DataRecord Ru_p3(
    R = Modelica.Constants.R/Ru_p3.MM,
    G_ref = 173385,
    H_ref = 135018,
    S_ref = -296.227,
    a1 = -1.1421e-05,
    a2 = -6041.95,
    a3 = 0.00047743,
    a4 = -91291,
    c1 = 72.3347,
    c2 = -320687,
    w_ref = 1110100,
    z = 3,
    a_0 = 3.72,
    MM = 0) "Ru+3";

  //504
  constant DataRecord RuCl_p1(
    R = Modelica.Constants.R/RuCl_p1.MM,
    G_ref = 21799,
    H_ref = -13891,
    S_ref = -45.187,
    a1 = 8.0785e-06,
    a2 = -1284.45,
    a3 = 0.00029131,
    a4 = -110960,
    c1 = 89.0962,
    c2 = -3301,
    w_ref = 300870,
    z = 1,
    a_0 = 3.72,
    MM = 0.13652) "RuCl+";

  //505
  constant DataRecord RuCl_p2(
    R = Modelica.Constants.R/RuCl_p2.MM,
    G_ref = 29706,
    H_ref = -33054,
    S_ref = -201.25,
    a1 = 8.08e-08,
    a2 = -3233.27,
    a3 = 0.00036706,
    a4 = -102901,
    c1 = 163.4768,
    c2 = 112391,
    w_ref = 746890,
    z = 2,
    a_0 = 3.72,
    MM = 0.13652) "RuCl+2";

  //506
  constant DataRecord RuCl2_aq(
    R = Modelica.Constants.R/RuCl2_aq.MM,
    G_ref = -104851,
    H_ref = -185770,
    S_ref = -20.083,
    a1 = 2.2288e-05,
    a2 = 2185.47,
    a3 = 0.00015485,
    a4 = -125302,
    c1 = 86.1138,
    c2 = 86680,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17197) "RuCl2,aq";

  //507
  constant DataRecord RuCl2_p1(
    R = Modelica.Constants.R/RuCl2_p1.MM,
    G_ref = -110625,
    H_ref = -213133,
    S_ref = -157.737,
    a1 = 1.3326e-05,
    a2 = -1.46,
    a3 = 0.00024055,
    a4 = -116265,
    c1 = 200.8797,
    c2 = 331293,
    w_ref = 469280,
    z = 1,
    a_0 = 3.72,
    MM = 0.17197) "RuCl2+";

  //508
  constant DataRecord RuCl3_aq(
    R = Modelica.Constants.R/RuCl3_aq.MM,
    G_ref = -245015,
    H_ref = -403254,
    S_ref = -168.197,
    a1 = 2.7501e-05,
    a2 = 3459.25,
    a3 = 0.00010459,
    a4 = -130570,
    c1 = 157.9251,
    c2 = 336276,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.20742) "RuCl3,aq";

  //509
  constant DataRecord RuCl3_n1(
    R = Modelica.Constants.R/RuCl3_n1.MM,
    G_ref = -227526,
    H_ref = -367481,
    S_ref = -41.422,
    a1 = 4.1857e-05,
    a2 = 6964.02,
    a3 = -3.3041e-05,
    a4 = -145059,
    c1 = 150.9089,
    c2 = 69701,
    w_ref = 743790,
    z = -1,
    a_0 = 3.72,
    MM = 0.20742) "RuCl3-";

  //510
  constant DataRecord RuCl4_n1(
    R = Modelica.Constants.R/RuCl4_n1.MM,
    G_ref = -375430,
    H_ref = -609065,
    S_ref = -244.346,
    a1 = 4.8708e-05,
    a2 = 8636.28,
    a3 = -9.865e-05,
    a4 = -151971,
    c1 = 195.9208,
    c2 = 127340,
    w_ref = 1052360,
    z = -1,
    a_0 = 3.72,
    MM = 0.24287) "RuCl4-";

  //511
  constant DataRecord RuCl4_n2(
    R = Modelica.Constants.R/RuCl4_n2.MM,
    G_ref = -351038,
    H_ref = -568187,
    S_ref = -123.01,
    a1 = 6.3475e-05,
    a2 = 12241.97,
    a3 = -0.00024031,
    a4 = -166879,
    c1 = 187.8599,
    c2 = -54241,
    w_ref = 1531970,
    z = -2,
    a_0 = 3.72,
    MM = 0.24287) "RuCl4-2";

  //512
  constant DataRecord RuCl5_n2(
    R = Modelica.Constants.R/RuCl5_n2.MM,
    G_ref = -505009,
    H_ref = -865900,
    S_ref = -494.047,
    a1 = 6.9051e-05,
    a2 = 13603.61,
    a3 = -0.0002939,
    a4 = -172506,
    c1 = 275.816,
    c2 = 71187,
    w_ref = 2094930,
    z = -2,
    a_0 = 3.72,
    MM = 0.27832) "RuCl5-2";

  //513
  constant DataRecord RuCl6_n3(
    R = Modelica.Constants.R/RuCl6_n3.MM,
    G_ref = -634043,
    H_ref = -1133048,
    S_ref = -780.609,
    a1 = 9.5454e-05,
    a2 = 20050.86,
    a3 = -0.00054736,
    a4 = -199158,
    c1 = 215.2484,
    c2 = -492252,
    w_ref = 3197120,
    z = -3,
    a_0 = 3.72,
    MM = 0.31377) "RuCl6-3";

  //514
  constant DataRecord RuO_aq(
    R = Modelica.Constants.R/RuO_aq.MM,
    G_ref = 962,
    H_ref = -59664,
    S_ref = -72.383,
    a1 = -3.9936e-06,
    a2 = -4228.31,
    a3 = 0.00040617,
    a4 = -98788,
    c1 = 34.9431,
    c2 = -91169,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.11707) "RuO,aq";

  //515
  constant DataRecord RuO_p1(
    R = Modelica.Constants.R/RuO_p1.MM,
    G_ref = -43723,
    H_ref = -96441,
    S_ref = -73.22,
    a1 = 9.7496e-06,
    a2 = -876.34,
    a3 = 0.00027527,
    a4 = -112646,
    c1 = -42.7977,
    c2 = -475545,
    w_ref = 344010,
    z = 1,
    a_0 = 3.72,
    MM = 0.11707) "RuO+";

  //516
  constant DataRecord S2_n2(
    R = Modelica.Constants.R/S2_n2.MM,
    G_ref = 79496,
    H_ref = 30125,
    S_ref = 28.451,
    a1 = 2.3345e-05,
    a2 = 2444.54,
    a3 = 0.0001445,
    a4 = -126378,
    c1 = -14.0147,
    c2 = -681804,
    w_ref = 1300510,
    z = -2,
    a_0 = 3.72,
    MM = 0.06412) "S2-2";

  //517
  constant DataRecord S2O3_n2(
    R = Modelica.Constants.R/S2O3_n2.MM,
    G_ref = -522582,
    H_ref = -648520,
    S_ref = 66.944,
    a1 = 2.7901e-05,
    a2 = 5227.95,
    a3 = -0.00032334,
    a4 = -137884,
    c1 = -0.2414,
    c2 = -615324,
    w_ref = 1242400,
    z = -2,
    a_0 = 4,
    MM = 0.11212) "S2O3-2";

  //518
  constant DataRecord HS2O3_n1(
    R = Modelica.Constants.R/HS2O3_n1.MM,
    G_ref = -532205,
    H_ref = -643918,
    S_ref = 127.612,
    a1 = 2.5926e-05,
    a2 = 3075.66,
    a3 = 0.00011945,
    a4 = -128984,
    c1 = 79.3893,
    c2 = -97136,
    w_ref = 488520,
    z = -1,
    a_0 = 3.72,
    MM = 0.11313) "HS2O3-";

  //519
  constant DataRecord H2S2O3_aq(
    R = Modelica.Constants.R/H2S2O3_aq.MM,
    G_ref = -535552,
    H_ref = -629274,
    S_ref = 188.28,
    a1 = 3.1063e-05,
    a2 = 4330.15,
    a3 = 7.0132e-05,
    a4 = -134168,
    c1 = 80.0705,
    c2 = 107412,
    w_ref = -142880,
    z = 0,
    a_0 = 3.72,
    MM = 0.11414) "H2S2O3,aq";

  //520
  constant DataRecord S2O4_n2(
    R = Modelica.Constants.R/S2O4_n2.MM,
    G_ref = -600404,
    H_ref = -753538,
    S_ref = 92.048,
    a1 = 2.7942e-05,
    a2 = 3568.12,
    a3 = 0.00010007,
    a4 = -131018,
    c1 = 13.8553,
    c2 = -553957,
    w_ref = 1203820,
    z = -2,
    a_0 = 5,
    MM = 0.12812) "S2O4-2";

  //521
  constant DataRecord HS2O4_n1(
    R = Modelica.Constants.R/HS2O4_n1.MM,
    G_ref = -614630,
    H_ref = -749354,
    S_ref = 152.716,
    a1 = 3.1865e-05,
    a2 = 4524.66,
    a3 = 6.2781e-05,
    a4 = -134976,
    c1 = 100.134,
    c2 = -12757,
    w_ref = 450200,
    z = -1,
    a_0 = 3.72,
    MM = 0.12913) "HS2O4-";

  //522
  constant DataRecord H2S2O4_aq(
    R = Modelica.Constants.R/H2S2O4_aq.MM,
    G_ref = -616722,
    H_ref = -733455,
    S_ref = 213.384,
    a1 = 3.6202e-05,
    a2 = 5584.05,
    a3 = 2.1016e-05,
    a4 = -139352,
    c1 = 100.844,
    c2 = 191786,
    w_ref = -180870,
    z = 0,
    a_0 = 3.72,
    MM = 0.13014) "H2S2O4,aq";

  //523
  constant DataRecord S2O5_n2(
    R = Modelica.Constants.R/S2O5_n2.MM,
    G_ref = -790776,
    H_ref = -970688,
    S_ref = 104.6,
    a1 = 3.0802e-05,
    a2 = 4265.38,
    a3 = 7.286e-05,
    a4 = -133901,
    c1 = 16.6168,
    c2 = -538615,
    w_ref = 1185870,
    z = -2,
    a_0 = 3.72,
    MM = 0.14412) "S2O5-2";

  //524
  constant DataRecord S2O6_n2(
    R = Modelica.Constants.R/S2O6_n2.MM,
    G_ref = -966504,
    H_ref = -1173194,
    S_ref = 125.52,
    a1 = 3.4416e-05,
    a2 = 5148.58,
    a3 = 3.802e-05,
    a4 = -137553,
    c1 = 21.0585,
    c2 = -513046,
    w_ref = 1154240,
    z = -2,
    a_0 = 4,
    MM = 0.16012) "S2O6-2";

  //525
  constant DataRecord S2O8_n2(
    R = Modelica.Constants.R/S2O8_n2.MM,
    G_ref = -1115036,
    H_ref = -1344738,
    S_ref = 244.346,
    a1 = 5.5907e-05,
    a2 = 10395.32,
    a3 = -0.000168,
    a4 = -159247,
    c1 = 54.238,
    c2 = -340038,
    w_ref = 974080,
    z = -2,
    a_0 = 4,
    MM = 0.19212) "S2O8-2";

  //526
  constant DataRecord S3_n2(
    R = Modelica.Constants.R/S3_n2.MM,
    G_ref = 73638,
    H_ref = 25941,
    S_ref = 66.107,
    a1 = 2.8309e-05,
    a2 = 3656.65,
    a3 = 9.686e-05,
    a4 = -131390,
    c1 = -1.5041,
    c2 = -620437,
    w_ref = 1244700,
    z = -2,
    a_0 = 3.72,
    MM = 0.09618) "S3-2";

  //527
  constant DataRecord S3O6_n2(
    R = Modelica.Constants.R/S3O6_n2.MM,
    G_ref = -958136,
    H_ref = -1167336,
    S_ref = 138.072,
    a1 = 3.5211e-05,
    a2 = 5342.59,
    a3 = 3.0409e-05,
    a4 = -138357,
    c1 = 23.7162,
    c2 = -497704,
    w_ref = 1135160,
    z = -2,
    a_0 = 3.72,
    MM = 0.19218) "S3O6-2";

  //528
  constant DataRecord S4_n2(
    R = Modelica.Constants.R/S4_n2.MM,
    G_ref = 69036,
    H_ref = 23012,
    S_ref = 103.345,
    a1 = 3.3213e-05,
    a2 = 4853.94,
    a3 = 4.9798e-05,
    a4 = -136340,
    c1 = 10.9123,
    c2 = -559074,
    w_ref = 1187840,
    z = -2,
    a_0 = 3.72,
    MM = 0.12824) "S4-2";

  //529
  constant DataRecord S4O6_n2(
    R = Modelica.Constants.R/S4O6_n2.MM,
    G_ref = -1040561,
    H_ref = -1224238,
    S_ref = 257.316,
    a1 = 4.2958e-05,
    a2 = 7234.22,
    a3 = -4.394e-05,
    a4 = -146176,
    c1 = 48.9704,
    c2 = -351966,
    w_ref = 954160,
    z = -2,
    a_0 = 3.72,
    MM = 0.22424) "S4O6-2";

  //530
  constant DataRecord S5_n2(
    R = Modelica.Constants.R/S5_n2.MM,
    G_ref = 65689,
    H_ref = 21338,
    S_ref = 140.582,
    a1 = 3.8119e-05,
    a2 = 6051.95,
    a3 = 2.715e-06,
    a4 = -141294,
    c1 = 23.163,
    c2 = -498561,
    w_ref = 1131810,
    z = -2,
    a_0 = 3.72,
    MM = 0.1603) "S5-2";

  //531
  constant DataRecord S5O6_n2(
    R = Modelica.Constants.R/S5O6_n2.MM,
    G_ref = -958136,
    H_ref = -1175704,
    S_ref = 167.36,
    a1 = 3.7123e-05,
    a2 = 5807.64,
    a3 = 1.2493e-05,
    a4 = -140277,
    c1 = 29.9474,
    c2 = -461909,
    w_ref = 1091020,
    z = -2,
    a_0 = 3.72,
    MM = 0.2563) "S5O6-2";

  //532
  constant DataRecord SEBACATE_aq(
    R = Modelica.Constants.R/SEBACATE_aq.MM,
    G_ref = -617726,
    H_ref = -1032946,
    S_ref = 249.785,
    a1 = 9.7178e-05,
    a2 = 19212.13,
    a3 = -0.00024418,
    a4 = -195694,
    c1 = 373.6877,
    c2 = -140758,
    w_ref = 966290,
    z = -2,
    a_0 = 3.72,
    MM = 0.20023) "SEBACATE,aq";

  //533
  constant DataRecord SEBACIC_ACID_aq(
    R = Modelica.Constants.R/SEBACIC_ACID_aq.MM,
    G_ref = -674628,
    H_ref = -1029264,
    S_ref = 453.127,
    a1 = 0.00010925,
    a2 = 21947.51,
    a3 = -0.00030604,
    a4 = -206999,
    c1 = 598.9701,
    c2 = 474507,
    w_ref = 54640,
    z = 0,
    a_0 = 3.72,
    MM = 0.20224) "SEBACIC-ACID,aq";

  //534
  constant DataRecord SERINE_aq(
    R = Modelica.Constants.R/SERINE_aq.MM,
    G_ref = -510866,
    H_ref = -714627,
    S_ref = 194.556,
    a1 = 4.1577e-05,
    a2 = 4914.28,
    a3 = 0.00047255,
    a4 = -136587,
    c1 = 121.8485,
    c2 = -90136,
    w_ref = -152380,
    z = 0,
    a_0 = 3.72,
    MM = 0.10509) "SERINE,aq";

  //535
  constant DataRecord SO2_aq(
    R = Modelica.Constants.R/SO2_aq.MM,
    G_ref = -301164,
    H_ref = -322980,
    S_ref = 161.921,
    a1 = 2.908e-05,
    a2 = 3844.68,
    a3 = 8.9466e-05,
    a4 = -132168,
    c1 = 130.5831,
    c2 = 270194,
    w_ref = -102970,
    z = 0,
    a_0 = 3.72,
    MM = 0.06406) "SO2,aq";

  //536
  constant DataRecord SO3_n2(
    R = Modelica.Constants.R/SO3_n2.MM,
    G_ref = -486599,
    H_ref = -635550,
    S_ref = -29.288,
    a1 = 1.0306e-05,
    a2 = -740.19,
    a3 = 0.00026984,
    a4 = -113211,
    c1 = -11.7014,
    c2 = -702255,
    w_ref = 1389510,
    z = -2,
    a_0 = 4.5,
    MM = 0.08006) "SO3-2";

  //537
  constant DataRecord SO4_n2(
    R = Modelica.Constants.R/SO4_n2.MM,
    G_ref = -744459,
    H_ref = -909602,
    S_ref = 18.828,
    a1 = 3.4733e-05,
    a2 = -830.36,
    a3 = -0.00025992,
    a4 = -112842,
    c1 = 6.8618,
    c2 = -753036,
    w_ref = 1316410,
    z = -2,
    a_0 = 4,
    MM = 0.09606) "SO4-2";

  //538
  constant DataRecord SUBERATE_aq(
    R = Modelica.Constants.R/SUBERATE_aq.MM,
    G_ref = -648394,
    H_ref = -999056,
    S_ref = 193.719,
    a1 = 7.9544e-05,
    a2 = 15215.49,
    a3 = -0.00015342,
    a4 = -179171,
    c1 = 208.4506,
    c2 = -162256,
    w_ref = 1050940,
    z = -2,
    a_0 = 3.72,
    MM = 0.17217) "SUBERATE,aq";

  //539
  constant DataRecord SUBERIC_ACID_aq(
    R = Modelica.Constants.R/SUBERIC_ACID_aq.MM,
    G_ref = -704962,
    H_ref = -994788,
    S_ref = 397.898,
    a1 = 9.1036e-05,
    a2 = 17818.48,
    a3 = -0.00021214,
    a4 = -189933,
    c1 = 453.8117,
    c2 = 299235,
    w_ref = 18030,
    z = 0,
    a_0 = 3.72,
    MM = 0.17419) "SUBERIC-ACID,aq";

  //540
  constant DataRecord SUCCINATE_aq(
    R = Modelica.Constants.R/SUCCINATE_aq.MM,
    G_ref = -687766,
    H_ref = -909392,
    S_ref = 81.588,
    a1 = 4.3755e-05,
    a2 = 6636.74,
    a3 = 0.00014944,
    a4 = -143704,
    c1 = -58.7362,
    c2 = -197401,
    w_ref = 1220470,
    z = -2,
    a_0 = 3.72,
    MM = 0.11607) "SUCCINATE,aq";

  //541
  constant DataRecord SUCCINIC_ACID_aq(
    R = Modelica.Constants.R/SUCCINIC_ACID_aq.MM,
    G_ref = -743915,
    H_ref = -912112,
    S_ref = 260.663,
    a1 = 5.4338e-05,
    a2 = 9502.24,
    a3 = -2.3527e-05,
    a4 = -155553,
    c1 = 213.6936,
    c2 = 12761,
    w_ref = -72970,
    z = 0,
    a_0 = 3.72,
    MM = 0.11809) "SUCCINIC-ACID,aq";

  //542
  constant DataRecord Sc_p3(
    R = Modelica.Constants.R/Sc_p3.MM,
    G_ref = -586597,
    H_ref = -614211,
    S_ref = -255.224,
    a1 = -8.832e-06,
    a2 = -5409.66,
    a3 = 0.00045258,
    a4 = -93906,
    c1 = 88.3385,
    c2 = -244580,
    w_ref = 1046130,
    z = 3,
    a_0 = 9,
    MM = 0.04496) "Sc+3";

  //543
  constant DataRecord ScOH_p2(
    R = Modelica.Constants.R/ScOH_p2.MM,
    G_ref = -79914,
    H_ref = -839729,
    S_ref = -65.689,
    a1 = -1.72e-06,
    a2 = -3672.63,
    a3 = 0.0003842,
    a4 = -101085,
    c1 = 28.37,
    c2 = -291453,
    w_ref = 541580,
    z = 2,
    a_0 = 3.72,
    MM = 0.061968) "ScOH+2";

  //544
  constant DataRecord ScO_p1(
    R = Modelica.Constants.R/ScO_p1.MM,
    G_ref = -768182,
    H_ref = -794123,
    S_ref = -15.062,
    a1 = 4.556e-07,
    a2 = -3140.89,
    a3 = 0.00036324,
    a4 = -103286,
    c1 = -89.861,
    c2 = -610207,
    w_ref = 253680,
    z = 1,
    a_0 = 3.72,
    MM = 0.06096) "ScO+";

  //545
  constant DataRecord HScO2_aq(
    R = Modelica.Constants.R/HScO2_aq.MM,
    G_ref = -969014,
    H_ref = -1022151,
    S_ref = 126.775,
    a1 = 1.3984e-05,
    a2 = 157.74,
    a3 = 0.00023458,
    a4 = -116922,
    c1 = -200.7019,
    c2 = -910208,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.077968) "HScO2,aq";

  //546
  constant DataRecord ScO2_n1(
    R = Modelica.Constants.R/ScO2_n1.MM,
    G_ref = -912530,
    H_ref = -979474,
    S_ref = 80.333,
    a1 = 1.2481e-05,
    a2 = -208.32,
    a3 = 0.0002488,
    a4 = -115407,
    c1 = -104.3653,
    c2 = -758505,
    w_ref = 559360,
    z = -1,
    a_0 = 3.72,
    MM = 0.07696) "ScO2-";

  //547
  constant DataRecord SeO3_n2(
    R = Modelica.Constants.R/SeO3_n2.MM,
    G_ref = -369866,
    H_ref = -509193,
    S_ref = 12.97,
    a1 = 1.5125e-05,
    a2 = 438.4,
    a3 = 0.0002231,
    a4 = -118081,
    c1 = -2.7267,
    c2 = -650265,
    w_ref = 1324570,
    z = -2,
    a_0 = 3.72,
    MM = 0.12696) "SeO3-2";

  //548
  constant DataRecord SeO4_n2(
    R = Modelica.Constants.R/SeO4_n2.MM,
    G_ref = -441412,
    H_ref = -599149,
    S_ref = 53.974,
    a1 = 2.4078e-05,
    a2 = 2622.45,
    a3 = 0.00013772,
    a4 = -127110,
    c1 = 5.8455,
    c2 = -600831,
    w_ref = 1263230,
    z = -2,
    a_0 = 4,
    MM = 0.14296) "SeO4-2";

  //549
  constant DataRecord SiF6_n2(
    R = Modelica.Constants.R/SiF6_n2.MM,
    G_ref = -2199529,
    H_ref = -2389064,
    S_ref = 122.173,
    a1 = 3.5694e-05,
    a2 = 5459.79,
    a3 = 2.5987e-05,
    a4 = -138846,
    c1 = 17.1418,
    c2 = -528393,
    w_ref = 1159640,
    z = -2,
    a_0 = 3.72,
    MM = 0.14208) "SiF6-2";

  //550
  constant DataRecord H4SiO4_aq(
    R = Modelica.Constants.R/H4SiO4_aq.MM,
    G_ref = -1309257,
    H_ref = -1458861,
    S_ref = 188.7,
    a1 = 7.836e-05,
    a2 = -8895.2,
    a3 = 0.00077906,
    a4 = -50210,
    c1 = 242.8,
    c2 = -869850,
    w_ref = 36360,
    z = 0,
    a_0 = 3.72,
    MM = 0.096112) "H4SiO4,aq";

  //551
  constant DataRecord SiO2_aq(
    R = Modelica.Constants.R/SiO2_aq.MM,
    G_ref = -834875,
    H_ref = -887757,
    S_ref = 46.56,
    a1 = 7.9496e-06,
    a2 = 711.28,
    a3 = 0.0008368,
    a4 = -112968,
    c1 = 134.8085,
    c2 = -1058134,
    w_ref = 143390,
    z = 0,
    a_0 = 3.72,
    MM = 0.06008) "SiO2,aq";

  //552
  constant DataRecord Sm_p2(
    R = Modelica.Constants.R/Sm_p2.MM,
    G_ref = -514632,
    H_ref = -504172,
    S_ref = -25.941,
    a1 = -1.477e-07,
    a2 = -3288.29,
    a3 = 0.000369,
    a4 = -102675,
    c1 = 43.9358,
    c2 = -218158,
    w_ref = 481660,
    z = 2,
    a_0 = 3.72,
    MM = 0.15035) "Sm+2";

  //553
  constant DataRecord Sm_p3(
    R = Modelica.Constants.R/Sm_p3.MM,
    G_ref = -665674,
    H_ref = -691197,
    S_ref = -212.129,
    a1 = -1.3416e-05,
    a2 = -6531.56,
    a3 = 0.0004973,
    a4 = -89274,
    c1 = 8.1107,
    c2 = -496005,
    w_ref = 960440,
    z = 3,
    a_0 = 9,
    MM = 0.15035) "Sm+3";

  //554
  constant DataRecord Sm_p4(
    R = Modelica.Constants.R/Sm_p4.MM,
    G_ref = -166942,
    H_ref = -235978,
    S_ref = -423.002,
    a1 = -1.7943e-05,
    a2 = -7635.8,
    a3 = 0.00054038,
    a4 = -84701,
    c1 = 168.2897,
    c2 = -128666,
    w_ref = 1551970,
    z = 4,
    a_0 = 3.72,
    MM = 0.15035) "Sm+4";

  //555
  constant DataRecord Sn_p2(
    R = Modelica.Constants.R/Sn_p2.MM,
    G_ref = -27489,
    H_ref = -8786,
    S_ref = -16.736,
    a1 = 3.93e-08,
    a2 = -3243.27,
    a3 = 0.0003674,
    a4 = -102859,
    c1 = 31.0285,
    c2 = -259069,
    w_ref = 469280,
    z = 2,
    a_0 = 6,
    MM = 0.11869) "Sn+2";

  //556
  constant DataRecord SnO_aq(
    R = Modelica.Constants.R/SnO_aq.MM,
    G_ref = -224262,
    H_ref = -251458,
    S_ref = 61.923,
    a1 = 1.1465e-05,
    a2 = -456.27,
    a3 = 0.00025847,
    a4 = -114382,
    c1 = -54.8025,
    c2 = -403103,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13469) "SnO,aq";

  //557
  constant DataRecord SnOH_p1(
    R = Modelica.Constants.R/SnOH_p1.MM,
    G_ref = -245182,
    H_ref = -266939,
    S_ref = 80.333,
    a1 = 1.2275e-05,
    a2 = -258.61,
    a3 = 0.00025074,
    a4 = -115198,
    c1 = 4.171,
    c2 = -236906,
    w_ref = 108570,
    z = 1,
    a_0 = 3.72,
    MM = 0.1357) "SnOH+";

  //558
  constant DataRecord HSnO2_n1(
    R = Modelica.Constants.R/HSnO2_n1.MM,
    G_ref = -407103,
    H_ref = -510866,
    S_ref = 39.33,
    a1 = 1.7561e-05,
    a2 = 1032.53,
    a3 = 0.00019991,
    a4 = -120537,
    c1 = -28.1884,
    c2 = -513900,
    w_ref = 622330,
    z = -1,
    a_0 = 3.72,
    MM = 0.1517) "HSnO2-";

  //559
  constant DataRecord Sr_CO3_aq(
    R = Modelica.Constants.R/Sr_CO3_aq.MM,
    G_ref = -1108174,
    H_ref = -1207586,
    S_ref = 35.564,
    a1 = -6.234e-07,
    a2 = -3407.45,
    a3 = 0.00037451,
    a4 = -102173,
    c1 = -55.4798,
    c2 = -404384,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.14763) "Sr(CO3),aq";

  //560
  constant DataRecord Sr_HCO3_p1(
    R = Modelica.Constants.R/Sr_HCO3_p1.MM,
    G_ref = -1157796,
    H_ref = -1220423,
    S_ref = 158.992,
    a1 = 1.5456e-05,
    a2 = 517.56,
    a3 = 0.00022025,
    a4 = -118407,
    c1 = 147.8835,
    c2 = 300453,
    w_ref = -9620,
    z = 1,
    a_0 = 3.72,
    MM = 0.14864) "Sr(HCO3)+";

  //561
  constant DataRecord Sr_p2(
    R = Modelica.Constants.R/Sr_p2.MM,
    G_ref = -563836,
    H_ref = -550907,
    S_ref = -31.506,
    a1 = 2.9585e-06,
    a2 = -4247.09,
    a3 = 0.00029299,
    a4 = -98717,
    c1 = 44.9579,
    c2 = -212623,
    w_ref = 475430,
    z = 2,
    a_0 = 5,
    MM = 0.08762) "Sr+2";

  //562
  constant DataRecord SrCl_p1(
    R = Modelica.Constants.R/SrCl_p1.MM,
    G_ref = -693707,
    H_ref = -710401,
    S_ref = 46.024,
    a1 = 1.1598e-05,
    a2 = -422.75,
    a3 = 0.00025691,
    a4 = -114524,
    c1 = 69.6226,
    c2 = -26054,
    w_ref = 160540,
    z = 1,
    a_0 = 3.72,
    MM = 0.12307) "SrCl+";

  //563
  constant DataRecord SrF_p1(
    R = Modelica.Constants.R/SrF_p1.MM,
    G_ref = -846381,
    H_ref = -881443,
    S_ref = -25.941,
    a1 = 9.774e-07,
    a2 = -3015.87,
    a3 = 0.00035882,
    a4 = -103805,
    c1 = 90.2514,
    c2 = 10339,
    w_ref = 270790,
    z = 1,
    a_0 = 3.72,
    MM = 0.10662) "SrF+";

  //564
  constant DataRecord SrOH_p1(
    R = Modelica.Constants.R/SrOH_p1.MM,
    G_ref = -725087,
    H_ref = -753120,
    S_ref = 61.086,
    a1 = 1.1975e-05,
    a2 = -331.46,
    a3 = 0.00025349,
    a4 = -114901,
    c1 = 19.9058,
    c2 = -191736,
    w_ref = 138320,
    z = 1,
    a_0 = 3.72,
    MM = 0.10463) "SrOH+";

  //565
  constant DataRecord THREONINE_aq(
    R = Modelica.Constants.R/THREONINE_aq.MM,
    G_ref = -513251,
    H_ref = -749354,
    S_ref = 222.589,
    a1 = 5.0731e-05,
    a2 = 6744.36,
    a3 = 0.00048757,
    a4 = -144151,
    c1 = 197.7024,
    c2 = -27631,
    w_ref = -194850,
    z = 0,
    a_0 = 3.72,
    MM = 0.11912) "THREONINE,aq";

  //566
  constant DataRecord TOLUENE_aq(
    R = Modelica.Constants.R/TOLUENE_aq.MM,
    G_ref = 126608,
    H_ref = 13724,
    S_ref = 183.678,
    a1 = 6.2868e-05,
    a2 = 9169.36,
    a3 = 0.00050772,
    a4 = -154176,
    c1 = 392.978,
    c2 = 121139,
    w_ref = -135900,
    z = 0,
    a_0 = 3.72,
    MM = 0.092133) "TOLUENE,aq";

  //567
  constant DataRecord TRYPTOPHAN_aq(
    R = Modelica.Constants.R/TRYPTOPHAN_aq.MM,
    G_ref = -112550,
    H_ref = -409195,
    S_ref = 248.948,
    a1 = 8.9527e-05,
    a2 = 14496.77,
    a3 = 0.00055192,
    a4 = -176201,
    c1 = 388.5777,
    c2 = 114349,
    w_ref = -89660,
    z = 0,
    a_0 = 3.72,
    MM = 0.20422) "TRYPTOPHAN,aq";

  //568
  constant DataRecord TYROSINE_aq(
    R = Modelica.Constants.R/TYROSINE_aq.MM,
    G_ref = -365263,
    H_ref = -658562,
    S_ref = 190.372,
    a1 = 7.7314e-05,
    a2 = 12055.23,
    a3 = 0.00053185,
    a4 = -166105,
    c1 = 279.0749,
    c2 = 32610,
    w_ref = -146020,
    z = 0,
    a_0 = 3.72,
    MM = 0.18118) "TYROSINE,aq";

  //569
  constant DataRecord Tb_p2(
    R = Modelica.Constants.R/Tb_p2.MM,
    G_ref = -332210,
    H_ref = -318821,
    S_ref = -12.134,
    a1 = 1.23e-07,
    a2 = -3223.44,
    a3 = 0.00036678,
    a4 = -102943,
    c1 = 40.4928,
    c2 = -223271,
    w_ref = 460240,
    z = 2,
    a_0 = 3.72,
    MM = 0.15892) "Tb+2";

  //570
  constant DataRecord Tb_p3(
    R = Modelica.Constants.R/Tb_p3.MM,
    G_ref = -667348,
    H_ref = -698310,
    S_ref = -225.936,
    a1 = -1.2236e-05,
    a2 = -6240.94,
    a3 = 0.00048526,
    a4 = -90471,
    c1 = 30.0633,
    c2 = -433785,
    w_ref = 1004450,
    z = 3,
    a_0 = 3.72,
    MM = 0.15892) "Tb+3";

  //571
  constant DataRecord Tb_p4(
    R = Modelica.Constants.R/Tb_p4.MM,
    G_ref = -369029,
    H_ref = -443086,
    S_ref = -435.973,
    a1 = -1.806e-05,
    a2 = -7662.74,
    a3 = 0.00054112,
    a4 = -84592,
    c1 = 167.3445,
    c2 = -137189,
    w_ref = 1568330,
    z = 4,
    a_0 = 3.72,
    MM = 0.15892) "Tb+4";

  //572
  constant DataRecord Th_p4(
    R = Modelica.Constants.R/Th_p4.MM,
    G_ref = -705004,
    H_ref = -769019,
    S_ref = -422.584,
    a1 = -1.7943e-05,
    a2 = -7635.8,
    a3 = 0.00054038,
    a4 = -84701,
    c1 = 168.2897,
    c2 = -128666,
    w_ref = 1551970,
    z = 4,
    a_0 = 11,
    MM = 0.23204) "Th+4";

  //573
  constant DataRecord Tl_p1(
    R = Modelica.Constants.R/Tl_p1.MM,
    G_ref = -32384,
    H_ref = 5356,
    S_ref = 125.52,
    a1 = 1.7943e-05,
    a2 = 1125.62,
    a3 = 0.00019627,
    a4 = -120922,
    c1 = 19.2564,
    c2 = -162758,
    w_ref = 40750,
    z = 1,
    a_0 = 2.5,
    MM = 0.20437) "Tl+";

  //574
  constant DataRecord Tl_p3(
    R = Modelica.Constants.R/Tl_p3.MM,
    G_ref = 214639,
    H_ref = 196648,
    S_ref = -192.464,
    a1 = -4.8564e-06,
    a2 = -4438.3,
    a3 = 0.00041429,
    a4 = -97922,
    c1 = 113.009,
    c2 = -128666,
    w_ref = 951940,
    z = 3,
    a_0 = 3.72,
    MM = 0.20437) "Tl+3";

  //575
  constant DataRecord TlOH_aq(
    R = Modelica.Constants.R/TlOH_aq.MM,
    G_ref = -274261,
    H_ref = -278905,
    S_ref = 135.98,
    a1 = 2.1656e-05,
    a2 = 2032.55,
    a3 = 0.00016055,
    a4 = -124671,
    c1 = -44.2588,
    c2 = -366456,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22138) "TlOH,aq";

  //576
  constant DataRecord TlOH_p2(
    R = Modelica.Constants.R/TlOH_p2.MM,
    G_ref = -18828,
    H_ref = -78659,
    S_ref = -99.161,
    a1 = 1.0021e-05,
    a2 = -809.23,
    a3 = 0.00027245,
    a4 = -112922,
    c1 = 55.4957,
    c2 = -213894,
    w_ref = 593790,
    z = 2,
    a_0 = 3.72,
    MM = 0.22138) "TlOH+2";

  //577
  constant DataRecord TlO_p1(
    R = Modelica.Constants.R/TlO_p1.MM,
    G_ref = -13389,
    H_ref = -21548,
    S_ref = 74.057,
    a1 = 1.2194e-05,
    a2 = -279.2,
    a3 = 0.00025171,
    a4 = -115114,
    c1 = -161.9007,
    c2 = -817315,
    w_ref = 118490,
    z = 1,
    a_0 = 3.72,
    MM = 0.22037) "TlO+";

  //578
  constant DataRecord HTlO2_aq(
    R = Modelica.Constants.R/HTlO2_aq.MM,
    G_ref = -240998,
    H_ref = -274470,
    S_ref = 222.17,
    a1 = 2.1656e-05,
    a2 = 2032.55,
    a3 = 0.00016055,
    a4 = -124671,
    c1 = -331.1523,
    c2 = -1363624,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.23738) "HTlO2,aq";

  //579
  constant DataRecord TlO2_n1(
    R = Modelica.Constants.R/TlO2_n1.MM,
    G_ref = -174054,
    H_ref = -219660,
    S_ref = 182.004,
    a1 = 2.1696e-05,
    a2 = 2043.05,
    a3 = 0.00015999,
    a4 = -124717,
    c1 = -241.1147,
    c2 = -1184645,
    w_ref = 405810,
    z = -1,
    a_0 = 3.72,
    MM = 0.23637) "TlO2-";

  //580
  constant DataRecord Tm_p2(
    R = Modelica.Constants.R/Tm_p2.MM,
    G_ref = -450198,
    H_ref = -441830,
    S_ref = -28.033,
    a1 = -1.941e-07,
    a2 = -3299.63,
    a3 = 0.00036943,
    a4 = -102629,
    c1 = 44.4701,
    c2 = -217304,
    w_ref = 484800,
    z = 2,
    a_0 = 3.72,
    MM = 0.16893) "Tm+2";

  //581
  constant DataRecord Tm_p3(
    R = Modelica.Constants.R/Tm_p3.MM,
    G_ref = -669022,
    H_ref = -705004,
    S_ref = -243.09,
    a1 = -1.3793e-05,
    a2 = -6623.77,
    a3 = 0.00050093,
    a4 = -88889,
    c1 = 35.4912,
    c2 = -419300,
    w_ref = 1018090,
    z = 3,
    a_0 = 3.72,
    MM = 0.16893) "Tm+3";

  //582
  constant DataRecord Tm_p4(
    R = Modelica.Constants.R/Tm_p4.MM,
    G_ref = -125938,
    H_ref = -202087,
    S_ref = -441.83,
    a1 = -1.8137e-05,
    a2 = -7681.66,
    a3 = 0.00054187,
    a4 = -84513,
    c1 = 167.3797,
    c2 = -140599,
    w_ref = 1579330,
    z = 4,
    a_0 = 3.72,
    MM = 0.16893) "Tm+4";

  //583
  constant DataRecord UNDECANOATE_aq(
    R = Modelica.Constants.R/UNDECANOATE_aq.MM,
    G_ref = -294595,
    H_ref = -704460,
    S_ref = 329.699,
    a1 = 0.00011152,
    a2 = 22461.59,
    a3 = -0.00031798,
    a4 = -209125,
    c1 = 851.5687,
    c2 = -72400,
    w_ref = 182260,
    z = -1,
    a_0 = 3.72,
    MM = 0.18528) "UNDECANOATE,aq";

  //584
  constant DataRecord UNDECANOIC_ACID_aq(
    R = Modelica.Constants.R/UNDECANOIC_ACID_aq.MM,
    G_ref = -322712,
    H_ref = -702368,
    S_ref = 430.952,
    a1 = 0.00011905,
    a2 = 24169.34,
    a3 = -0.00035678,
    a4 = -216187,
    c1 = 818.3176,
    c2 = 747338,
    w_ref = 39960,
    z = 0,
    a_0 = 3.72,
    MM = 0.18628) "UNDECANOIC-ACID,aq";

  //585
  constant DataRecord UREA_aq(
    R = Modelica.Constants.R/UREA_aq.MM,
    G_ref = -203844,
    H_ref = -318402,
    S_ref = 176.983,
    a1 = 3.2283e-05,
    a2 = 3055.62,
    a3 = 0.00045753,
    a4 = -128901,
    c1 = 98.6755,
    c2 = -110219,
    w_ref = -125770,
    z = 0,
    a_0 = 3.72,
    MM = 0.060055) "UREA,aq";

  //586
  constant DataRecord U_p3(
    R = Modelica.Constants.R/U_p3.MM,
    G_ref = -475093,
    H_ref = -489110,
    S_ref = -192.882,
    a1 = -1.1898e-05,
    a2 = -6159.68,
    a3 = 0.00048233,
    a4 = -90805,
    c1 = 24.2442,
    c2 = -437195,
    w_ref = 951940,
    z = 3,
    a_0 = 3.72,
    MM = 0.23803) "U+3";

  //587
  constant DataRecord U_p4(
    R = Modelica.Constants.R/U_p4.MM,
    G_ref = -529904,
    H_ref = -591199,
    S_ref = -416.726,
    a1 = -1.7923e-05,
    a2 = -7628.23,
    a3 = 0.00053955,
    a4 = -84734,
    c1 = 168.2776,
    c2 = -125261,
    w_ref = 1541180,
    z = 4,
    a_0 = 3.72,
    MM = 0.23803) "U+4";

  //588
  constant DataRecord UO2_p1(
    R = Modelica.Constants.R/UO2_p1.MM,
    G_ref = -961023,
    H_ref = -1025080,
    S_ref = -25.104,
    a1 = 1.4128e-05,
    a2 = 193.05,
    a3 = 0.00023315,
    a4 = -117068,
    c1 = -3.7664,
    c2 = -315319,
    w_ref = 267270,
    z = 1,
    a_0 = 3.72,
    MM = 0.27003) "UO2+";

  //589
  constant DataRecord UO2_p2(
    R = Modelica.Constants.R/UO2_p2.MM,
    G_ref = -952613,
    H_ref = -1019013,
    S_ref = -98.324,
    a1 = 1.2659e-05,
    a2 = -1718.95,
    a3 = 0.00064152,
    a4 = -109165,
    c1 = 88.1569,
    c2 = 17991,
    w_ref = 589900,
    z = 2,
    a_0 = 3.72,
    MM = 0.27003) "UO2+2";

  //590
  constant DataRecord VALINE_aq(
    R = Modelica.Constants.R/VALINE_aq.MM,
    G_ref = -356895,
    H_ref = -616303,
    S_ref = 178.238,
    a1 = 5.8934e-05,
    a2 = 8382.35,
    a3 = 0.00050142,
    a4 = -150921,
    c1 = 283.2945,
    c2 = 34589,
    w_ref = -127650,
    z = 0,
    a_0 = 3.72,
    MM = 0.11714) "VALINE,aq";

  //591
  constant DataRecord VO_p1(
    R = Modelica.Constants.R/VO_p1.MM,
    G_ref = -443922,
    H_ref = -457730,
    S_ref = 20.502,
    a1 = -1.3005e-05,
    a2 = -6428.34,
    a3 = 0.0004925,
    a4 = -89692,
    c1 = -118.5202,
    c2 = -692879,
    w_ref = 200790,
    z = 1,
    a_0 = 3.72,
    MM = 0.06694) "VO+";

  //592
  constant DataRecord VO_p2(
    R = Modelica.Constants.R/VO_p2.MM,
    G_ref = -446433,
    H_ref = -486599,
    S_ref = -133.888,
    a1 = -5.4271e-06,
    a2 = -4578.59,
    a3 = 0.00041998,
    a4 = -97341,
    c1 = 162.4484,
    c2 = 140649,
    w_ref = 647470,
    z = 2,
    a_0 = 3.72,
    MM = 0.06694) "VO+2";

  //593
  constant DataRecord VOH_p1(
    R = Modelica.Constants.R/VOH_p1.MM,
    G_ref = -417563,
    H_ref = -477394,
    S_ref = -68.618,
    a1 = -6.3693e-06,
    a2 = -4809.51,
    a3 = 0.00042924,
    a4 = -96387,
    c1 = 124.5698,
    c2 = 109115,
    w_ref = 334800,
    z = 1,
    a_0 = 3.72,
    MM = 0.067948) "VOH+";

  //594
  constant DataRecord VOH_p2(
    R = Modelica.Constants.R/VOH_p2.MM,
    G_ref = -466516,
    H_ref = -499570,
    S_ref = -44.35,
    a1 = -5.075e-07,
    a2 = -3376.78,
    a3 = 0.00037262,
    a4 = -102311,
    c1 = 11.3014,
    c2 = -340887,
    w_ref = 510700,
    z = 2,
    a_0 = 3.72,
    MM = 0.067948) "VOH+2";

  //595
  constant DataRecord VO2_p1(
    R = Modelica.Constants.R/VO2_p1.MM,
    G_ref = -587015,
    H_ref = -649775,
    S_ref = -42.258,
    a1 = 1.3642e-05,
    a2 = 74.89,
    a3 = 0.00023772,
    a4 = -116579,
    c1 = 19.2012,
    c2 = -243726,
    w_ref = 293010,
    z = 1,
    a_0 = 3.72,
    MM = 0.08294) "VO2+";

  //596
  constant DataRecord VOOH_p1(
    R = Modelica.Constants.R/VOOH_p1.MM,
    G_ref = -651449,
    H_ref = -743078,
    S_ref = -74.057,
    a1 = -6.5672e-06,
    a2 = -4856.87,
    a3 = 0.00043089,
    a4 = -96190,
    c1 = 129.0944,
    c2 = 121897,
    w_ref = 344010,
    z = 1,
    a_0 = 3.72,
    MM = 0.083948) "VOOH+";

  //597
  constant DataRecord WO4_n2(
    R = Modelica.Constants.R/WO4_n2.MM,
    G_ref = -914204,
    H_ref = -989516,
    S_ref = 40.585,
    a1 = 3.0156e-05,
    a2 = 3679.12,
    a3 = 0.00018782,
    a4 = -131478,
    c1 = 34.8561,
    c2 = -506226,
    w_ref = 1282690,
    z = -2,
    a_0 = 5,
    MM = 0.24785) "WO4-2";

  //598
  constant DataRecord HWO4_n1(
    R = Modelica.Constants.R/HWO4_n1.MM,
    G_ref = -934706,
    H_ref = -1066920,
    S_ref = 130.541,
    a1 = 3.2493e-05,
    a2 = 4678.55,
    a3 = 5.6618e-05,
    a4 = -135612,
    c1 = 145.3911,
    c2 = 133829,
    w_ref = 483630,
    z = -1,
    a_0 = 3.72,
    MM = 0.24886) "HWO4-";

  //599
  constant DataRecord Xe_aq(
    R = Modelica.Constants.R/Xe_aq.MM,
    G_ref = 13493,
    H_ref = -18870,
    S_ref = 61.17,
    a1 = 3.1462e-05,
    a2 = 4426.38,
    a3 = 6.6605e-05,
    a4 = -134570,
    c1 = 176.1615,
    c2 = 425312,
    w_ref = -92630,
    z = 0,
    a_0 = 3.72,
    MM = 0.1313) "Xe,aq";

  //600
  constant DataRecord Y_p3(
    R = Modelica.Constants.R/Y_p3.MM,
    G_ref = -685339,
    H_ref = -715046,
    S_ref = -251.04,
    a1 = -8.5617e-06,
    a2 = -5344.68,
    a3 = 0.00045024,
    a4 = -94173,
    c1 = 90.1091,
    c2 = -236906,
    w_ref = 1041400,
    z = 3,
    a_0 = 9,
    MM = 0.08891) "Y+3";

  //601
  constant DataRecord YOH_p2(
    R = Modelica.Constants.R/YOH_p2.MM,
    G_ref = -878640,
    H_ref = -924246,
    S_ref = -71.965,
    a1 = 1.0454e-05,
    a2 = -703.16,
    a3 = 0.00026819,
    a4 = -113361,
    c1 = 33.5289,
    c2 = -276964,
    w_ref = 552330,
    z = 2,
    a_0 = 3.72,
    MM = 0.10592) "YOH+2";

  //602
  constant DataRecord YO_p1(
    R = Modelica.Constants.R/YO_p1.MM,
    G_ref = -828850,
    H_ref = -856046,
    S_ref = -9.205,
    a1 = 1.0786e-05,
    a2 = -623,
    a3 = 0.00026524,
    a4 = -113692,
    c1 = -94.6793,
    c2 = -623843,
    w_ref = 243970,
    z = 1,
    a_0 = 3.72,
    MM = 0.10491) "YO+";

  //603
  constant DataRecord HYO2_aq(
    R = Modelica.Constants.R/HYO2_aq.MM,
    G_ref = -1011273,
    H_ref = -1065246,
    S_ref = 133.051,
    a1 = 1.8621e-05,
    a2 = 1289.84,
    a3 = 0.00019013,
    a4 = -121600,
    c1 = -211.9807,
    c2 = -949417,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.12192) "HYO2,aq";

  //604
  constant DataRecord YO2_n1(
    R = Modelica.Constants.R/YO2_n1.MM,
    G_ref = -951442,
    H_ref = -1019222,
    S_ref = 87.027,
    a1 = 1.892e-05,
    a2 = 1363.11,
    a3 = 0.00018715,
    a4 = -121905,
    c1 = -113.0751,
    c2 = -785776,
    w_ref = 549990,
    z = -1,
    a_0 = 3.72,
    MM = 0.12091) "YO2-";

  //605
  constant DataRecord Yb_p2(
    R = Modelica.Constants.R/Yb_p2.MM,
    G_ref = -537644,
    H_ref = -530531,
    S_ref = -46.861,
    a1 = -5.535e-07,
    a2 = -3388.2,
    a3 = 0.0003731,
    a4 = -102261,
    c1 = 49.1256,
    c2 = -210489,
    w_ref = 514000,
    z = 2,
    a_0 = 3.72,
    MM = 0.17304) "Yb+2";

  //606
  constant DataRecord Yb_p3(
    R = Modelica.Constants.R/Yb_p3.MM,
    G_ref = -640152,
    H_ref = -670695,
    S_ref = -238.07,
    a1 = -1.4637e-05,
    a2 = -6829.67,
    a3 = 0.00050902,
    a4 = -88040,
    c1 = 30.7662,
    c2 = -437199,
    w_ref = 1022700,
    z = 3,
    a_0 = 3.72,
    MM = 0.17304) "Yb+3";

  //607
  constant DataRecord Yb_p4(
    R = Modelica.Constants.R/Yb_p4.MM,
    G_ref = 15899,
    H_ref = -56902,
    S_ref = -446.014,
    a1 = -1.8176e-05,
    a2 = -7692.45,
    a3 = 0.00054259,
    a4 = -84467,
    c1 = 167.1554,
    c2 = -143156,
    w_ref = 1584900,
    z = 4,
    a_0 = 3.72,
    MM = 0.17304) "Yb+4";

  //608
  constant DataRecord Zn_p2(
    R = Modelica.Constants.R/Zn_p2.MM,
    G_ref = -147277,
    H_ref = -153385,
    S_ref = -109.621,
    a1 = -4.4668e-06,
    a2 = -3735.89,
    a3 = 0.0002564,
    a4 = -100826,
    c1 = 78.4082,
    c2 = -224681,
    w_ref = 609780,
    z = 2,
    a_0 = 6,
    MM = 0.06537) "Zn+2";

  //609
  constant DataRecord ZnCl_p1(
    R = Modelica.Constants.R/ZnCl_p1.MM,
    G_ref = -279700,
    H_ref = -277148,
    S_ref = 96.232,
    a1 = 6.9383e-06,
    a2 = -1560.34,
    a3 = 0.00030162,
    a4 = -109822,
    c1 = 82.4026,
    c2 = 42639,
    w_ref = 84730,
    z = 1,
    a_0 = 3.72,
    MM = 0.10082) "ZnCl+";

  //610
  constant DataRecord ZnCl2_aq(
    R = Modelica.Constants.R/ZnCl2_aq.MM,
    G_ref = -411287,
    H_ref = -456391,
    S_ref = 113.094,
    a1 = 2.1542e-05,
    a2 = 2005.35,
    a3 = 0.00016147,
    a4 = -124562,
    c1 = 109.4233,
    c2 = 168774,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.13627) "ZnCl2,aq";

  //611
  constant DataRecord ZnCl3_n1(
    R = Modelica.Constants.R/ZnCl3_n1.MM,
    G_ref = -541033,
    H_ref = -632035,
    S_ref = 104.6,
    a1 = 4.0014e-05,
    a2 = 6515.83,
    a3 = -1.5811e-05,
    a4 = -143206,
    c1 = 176.9464,
    c2 = 230735,
    w_ref = 523540,
    z = -1,
    a_0 = 3.72,
    MM = 0.17172) "ZnCl3-";

  //612
  constant DataRecord ZnOH_p1(
    R = Modelica.Constants.R/ZnOH_p1.MM,
    G_ref = -339699,
    H_ref = -363966,
    S_ref = 62.76,
    a1 = 4.8112e-06,
    a2 = -2078.49,
    a3 = 0.00032173,
    a4 = -107675,
    c1 = 62.888,
    c2 = -41735,
    w_ref = 136400,
    z = 1,
    a_0 = 3.72,
    MM = 0.082378) "ZnOH+";

  //613
  constant DataRecord ZnO_aq(
    R = Modelica.Constants.R/ZnO_aq.MM,
    G_ref = -282085,
    H_ref = -327565,
    S_ref = -8.368,
    a1 = -5.1957e-06,
    a2 = -4521.73,
    a3 = 0.00041764,
    a4 = -97575,
    c1 = 0.1234,
    c2 = -212192,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.08137) "ZnO,aq";

  //614
  constant DataRecord HZnO2_n1(
    R = Modelica.Constants.R/HZnO2_n1.MM,
    G_ref = -463252,
    H_ref = -595676,
    S_ref = -66.944,
    a1 = 2.3527e-06,
    a2 = -2679.73,
    a3 = 0.00034555,
    a4 = -105190,
    c1 = 146.8057,
    c2 = 43488,
    w_ref = 781110,
    z = -1,
    a_0 = 3.72,
    MM = 0.098378) "HZnO2-";

  //615
  constant DataRecord ZnO2_n2(
    R = Modelica.Constants.R/ZnO2_n2.MM,
    G_ref = -390325,
    H_ref = -552706,
    S_ref = -167.36,
    a1 = -2.3259e-06,
    a2 = -3821.96,
    a3 = 0.00039037,
    a4 = -100470,
    c1 = 136.3302,
    c2 = -254806,
    w_ref = 1598960,
    z = -2,
    a_0 = 3.72,
    MM = 0.09737) "ZnO2-2";

  //616
  constant DataRecord Bi_Ac_p2(
    R = Modelica.Constants.R/Bi_Ac_p2.MM,
    G_ref = -283466,
    H_ref = -403421,
    S_ref = -63.597,
    a1 = 2.1627e-05,
    a2 = 2024.47,
    a3 = 0.00016112,
    a4 = -124637,
    c1 = 451.6482,
    c2 = 1180875,
    w_ref = 538060,
    z = 2,
    a_0 = 3.72,
    MM = 0.26802) "Bi(Ac)+2";

  //617
  constant DataRecord Bi_Ac2_p1(
    R = Modelica.Constants.R/Bi_Ac2_p1.MM,
    G_ref = -657181,
    H_ref = -888598,
    S_ref = 40.166,
    a1 = 5.0108e-05,
    a2 = 8979.66,
    a3 = -0.00011245,
    a4 = -153390,
    c1 = 843.8919,
    c2 = 2662292,
    w_ref = 169280,
    z = 1,
    a_0 = 3.72,
    MM = 0.32707) "Bi(Ac)2+";

  //618
  constant DataRecord Bi_Ac3_aq(
    R = Modelica.Constants.R/Bi_Ac3_aq.MM,
    G_ref = -1027339,
    H_ref = -1377540,
    S_ref = 119.662,
    a1 = 8.2615e-05,
    a2 = 16916.29,
    a3 = -0.00042428,
    a4 = -186201,
    c1 = 1312.1693,
    c2 = 4348138,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.29011) "Bi(Ac)3,aq";

  //619
  constant DataRecord Dy_Ac_p2(
    R = Modelica.Constants.R/Dy_Ac_p2.MM,
    G_ref = -1048469,
    H_ref = -1197252,
    S_ref = -143.093,
    a1 = 1.1962e-05,
    a2 = -333.46,
    a3 = 0.00025339,
    a4 = -114893,
    c1 = 268.0772,
    c2 = 503570,
    w_ref = 660650,
    z = 2,
    a_0 = 3.72,
    MM = 0.22154) "Dy(Ac)+2";

  //620
  constant DataRecord Dy_Ac2_p1(
    R = Modelica.Constants.R/Dy_Ac2_p1.MM,
    G_ref = -1428627,
    H_ref = -1697491,
    S_ref = -68.199,
    a1 = 3.9437e-05,
    a2 = 6373.49,
    a3 = -9.887e-06,
    a4 = -142616,
    c1 = 478.9333,
    c2 = 1340792,
    w_ref = 334800,
    z = 1,
    a_0 = 3.72,
    MM = 0.28059) "Dy(Ac)2+";

  //621
  constant DataRecord Dy_Ac3_aq(
    R = Modelica.Constants.R/Dy_Ac3_aq.MM,
    G_ref = -1805773,
    H_ref = -2203930,
    S_ref = -22.594,
    a1 = 7.0102e-05,
    a2 = 13860.71,
    a3 = -0.00030408,
    a4 = -173569,
    c1 = 669.5417,
    c2 = 2114527,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.24363) "Dy(Ac)3,aq";

  //622
  constant DataRecord Ho_Ac_p2(
    R = Modelica.Constants.R/Ho_Ac_p2.MM,
    G_ref = -1059682,
    H_ref = -1207168,
    S_ref = -137.654,
    a1 = 1.1359e-05,
    a2 = -483.38,
    a3 = 0.00025978,
    a4 = -114269,
    c1 = 258.4369,
    c2 = 472888,
    w_ref = 651830,
    z = 2,
    a_0 = 3.72,
    MM = 0.22397) "Ho(Ac)+2";

  //623
  constant DataRecord Ho_Ac2_p1(
    R = Modelica.Constants.R/Ho_Ac2_p1.MM,
    G_ref = -1439798,
    H_ref = -1706779,
    S_ref = -61.086,
    a1 = 3.8753e-05,
    a2 = 6205.58,
    a3 = -3.113e-06,
    a4 = -141921,
    c1 = 460.4906,
    c2 = 1280927,
    w_ref = 321580,
    z = 1,
    a_0 = 3.72,
    MM = 0.28302) "Ho(Ac)2+";

  //624
  constant DataRecord Ho_Ac3_aq(
    R = Modelica.Constants.R/Ho_Ac3_aq.MM,
    G_ref = -1816902,
    H_ref = -2212231,
    S_ref = -12.97,
    a1 = 6.9389e-05,
    a2 = 13687.12,
    a3 = -0.00029736,
    a4 = -172854,
    c1 = 640.4302,
    c2 = 2013345,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.24606) "Ho(Ac)3,aq";

  //625
  constant DataRecord Er_Ac_p2(
    R = Modelica.Constants.R/Er_Ac_p2.MM,
    G_ref = -1053406,
    H_ref = -1207168,
    S_ref = -161.502,
    a1 = 1.0588e-05,
    a2 = -670.36,
    a3 = 0.00026689,
    a4 = -113499,
    c1 = 256.2512,
    c2 = 453709,
    w_ref = 688020,
    z = 2,
    a_0 = 3.72,
    MM = 0.2263) "Er(Ac)+2";

  //626
  constant DataRecord Er_Ac2_p1(
    R = Modelica.Constants.R/Er_Ac2_p1.MM,
    G_ref = -1433522,
    H_ref = -1709331,
    S_ref = -93.303,
    a1 = 3.7933e-05,
    a2 = 6005.21,
    a3 = 4.795e-06,
    a4 = -141093,
    c1 = 454.505,
    c2 = 1243510,
    w_ref = 373460,
    z = 1,
    a_0 = 3.72,
    MM = 0.28535) "Er(Ac)2+";

  //627
  constant DataRecord Er_Ac3_aq(
    R = Modelica.Constants.R/Er_Ac3_aq.MM,
    G_ref = -1810626,
    H_ref = -2217487,
    S_ref = -55.229,
    a1 = 6.8281e-05,
    a2 = 13416.29,
    a3 = -0.00028669,
    a4 = -171732,
    c1 = 622.2365,
    c2 = 1950104,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.24839) "Er(Ac)3,aq";

  //628
  constant DataRecord Tm_Ac_p2(
    R = Modelica.Constants.R/Tm_Ac_p2.MM,
    G_ref = -1053406,
    H_ref = -1207084,
    S_ref = -160.247,
    a1 = 1.0572e-05,
    a2 = -674.13,
    a3 = 0.000267,
    a4 = -113483,
    c1 = 255.8227,
    c2 = 453709,
    w_ref = 683330,
    z = 2,
    a_0 = 3.72,
    MM = 0.22797) "Tm(Ac)+2";

  //629
  constant DataRecord Tm_Ac2_p1(
    R = Modelica.Constants.R/Tm_Ac2_p1.MM,
    G_ref = -1433522,
    H_ref = -1709122,
    S_ref = -91.63,
    a1 = 3.7916e-05,
    a2 = 6001.86,
    a3 = 4.732e-06,
    a4 = -141080,
    c1 = 454.0335,
    c2 = 1243510,
    w_ref = 368320,
    z = 1,
    a_0 = 3.72,
    MM = 0.28702) "Tm(Ac)2+";

  //630
  constant DataRecord Tm_Ac3_aq(
    R = Modelica.Constants.R/Tm_Ac3_aq.MM,
    G_ref = -1810626,
    H_ref = -2217102,
    S_ref = -53.137,
    a1 = 6.8281e-05,
    a2 = 13416.29,
    a3 = -0.00028669,
    a4 = -171732,
    c1 = 622.2365,
    c2 = 1950104,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.25006) "Tm(Ac)3,aq";

  //631
  constant DataRecord Be_Ac_p1(
    R = Modelica.Constants.R/Be_Ac_p1.MM,
    G_ref = -728100,
    H_ref = -891359,
    S_ref = -190.79,
    a1 = 2.1253e-05,
    a2 = 1934.6,
    a3 = 0.00016428,
    a4 = -124265,
    c1 = 311.9046,
    c2 = 700443,
    w_ref = 521540,
    z = 1,
    a_0 = 3.72,
    MM = 0.068054) "Be(Ac)+";

  //632
  constant DataRecord Be_Ac2_aq(
    R = Modelica.Constants.R/Be_Ac2_aq.MM,
    G_ref = -1103488,
    H_ref = -1406786,
    S_ref = -183.259,
    a1 = 4.9138e-05,
    a2 = 8742.93,
    a3 = -0.00010318,
    a4 = -152411,
    c1 = 552.0181,
    c2 = 1706043,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.063097) "Be(Ac)2,aq";

  //633
  constant DataRecord Ra_Ac_p1(
    R = Modelica.Constants.R/Ra_Ac_p1.MM,
    G_ref = -936798,
    H_ref = -1001566,
    S_ref = 200.832,
    a1 = 2.8683e-05,
    a2 = 3747.44,
    a3 = 9.3366e-05,
    a4 = -131763,
    c1 = 184.8182,
    c2 = 449236,
    w_ref = -73390,
    z = 1,
    a_0 = 3.72,
    MM = 0.28504) "Ra(Ac)+";

  //634
  constant DataRecord Ra_Ac2_aq(
    R = Modelica.Constants.R/Ra_Ac2_aq.MM,
    G_ref = -1309341,
    H_ref = -1478040,
    S_ref = 329.699,
    a1 = 5.9657e-05,
    a2 = 11310.36,
    a3 = -0.00020387,
    a4 = -163025,
    c1 = 411.0006,
    c2 = 1215908,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.28009) "Ra(Ac)2,aq";

  //635
  constant DataRecord Au_Ac_aq(
    R = Modelica.Constants.R/Au_Ac_aq.MM,
    G_ref = -208656,
    H_ref = -285809,
    S_ref = 200.832,
    a1 = 4.2731e-05,
    a2 = 7178.74,
    a3 = -4.171e-05,
    a4 = -145946,
    c1 = 162.1015,
    c2 = 350799,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22401) "Au(Ac),aq";

  //636
  constant DataRecord Au_Ac2_n1(
    R = Modelica.Constants.R/Au_Ac2_n1.MM,
    G_ref = -578396,
    H_ref = -781362,
    S_ref = 256.479,
    a1 = 7.6114e-05,
    a2 = 15329.84,
    a3 = -0.00036206,
    a4 = -179644,
    c1 = 378.4888,
    c2 = 1004968,
    w_ref = 293300,
    z = -1,
    a_0 = 3.72,
    MM = 0.31506) "Au(Ac)2-";

  //637
  constant DataRecord Li_Ac_aq(
    R = Modelica.Constants.R/Li_Ac_aq.MM,
    G_ref = -663624,
    H_ref = -770860,
    S_ref = 81.588,
    a1 = 3.5095e-05,
    a2 = 5312.68,
    a3 = 3.1962e-05,
    a4 = -138231,
    c1 = 237.1353,
    c2 = 611596,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.033984) "Li(Ac),aq";

  //638
  constant DataRecord Li_Ac2_n1(
    R = Modelica.Constants.R/Li_Ac2_n1.MM,
    G_ref = -1032653,
    H_ref = -1274739,
    S_ref = 106.692,
    a1 = 6.8372e-05,
    a2 = 13439.47,
    a3 = -0.0002878,
    a4 = -171829,
    c1 = 545.7497,
    c2 = 1513813,
    w_ref = 519740,
    z = -1,
    a_0 = 3.72,
    MM = 0.12503) "Li(Ac)2-";

  //639
  constant DataRecord Na_Ac_aq(
    R = Modelica.Constants.R/Na_Ac_aq.MM,
    G_ref = -630612,
    H_ref = -726091,
    S_ref = 143.093,
    a1 = 3.4942e-05,
    a2 = 5277.07,
    a3 = 3.2987e-05,
    a4 = -138085,
    c1 = 208.777,
    c2 = 513030,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.050034) "Na(Ac),aq";

  //640
  constant DataRecord Na_Ac2_n1(
    R = Modelica.Constants.R/Na_Ac2_n1.MM,
    G_ref = -997758,
    H_ref = -1223402,
    S_ref = 184.096,
    a1 = 6.7807e-05,
    a2 = 13300.27,
    a3 = -0.00028207,
    a4 = -171251,
    c1 = 479.6692,
    c2 = 1321500,
    w_ref = 403040,
    z = -1,
    a_0 = 3.72,
    MM = 0.14108) "Na(Ac)2-";

  //641
  constant DataRecord K_Ac_aq(
    R = Modelica.Constants.R/K_Ac_aq.MM,
    G_ref = -650277,
    H_ref = -733120,
    S_ref = 199.158,
    a1 = 4.143e-05,
    a2 = 6859.12,
    a3 = -2.8761e-05,
    a4 = -144624,
    c1 = 169.716,
    c2 = 377259,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.066144) "K(Ac),aq";

  //642
  constant DataRecord K_Ac2_n1(
    R = Modelica.Constants.R/K_Ac2_n1.MM,
    G_ref = -1016670,
    H_ref = -1225494,
    S_ref = 253.969,
    a1 = 7.4676e-05,
    a2 = 14978.05,
    a3 = -0.00034808,
    a4 = -178188,
    c1 = 393.6793,
    c2 = 1056598,
    w_ref = 296940,
    z = -1,
    a_0 = 3.72,
    MM = 0.15719) "K(Ac)2-";

  //643
  constant DataRecord Mg_Ac_p1(
    R = Modelica.Constants.R/Mg_Ac_p1.MM,
    G_ref = -830608,
    H_ref = -960144,
    S_ref = -54.81,
    a1 = 2.3004e-05,
    a2 = 2360.78,
    a3 = 0.00014787,
    a4 = -126030,
    c1 = 270.4107,
    c2 = 622972,
    w_ref = 313090,
    z = 1,
    a_0 = 3.72,
    MM = 0.083354) "Mg(Ac)+";

  //644
  constant DataRecord Mg_Ac2_aq(
    R = Modelica.Constants.R/Mg_Ac2_aq.MM,
    G_ref = -1204281,
    H_ref = -1461584,
    S_ref = -5.021,
    a1 = 5.1874e-05,
    a2 = 9409.73,
    a3 = -0.00012909,
    a4 = -155168,
    c1 = 508.5288,
    c2 = 1554887,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.078397) "Mg(Ac)2,aq";

  //645
  constant DataRecord Sr_Ac_p1(
    R = Modelica.Constants.R/Sr_Ac_p1.MM,
    G_ref = -939350,
    H_ref = -1034368,
    S_ref = 83.68,
    a1 = 2.4937e-05,
    a2 = 2833.32,
    a3 = 0.00012922,
    a4 = -127980,
    c1 = 225.1448,
    c2 = 532652,
    w_ref = 103850,
    z = 1,
    a_0 = 3.72,
    MM = 0.14666) "Sr(Ac)+";

  //646
  constant DataRecord Sr_Ac2_aq(
    R = Modelica.Constants.R/Sr_Ac2_aq.MM,
    G_ref = -1312144,
    H_ref = -1521888,
    S_ref = 176.565,
    a1 = 5.4817e-05,
    a2 = 10129.55,
    a3 = -0.00015767,
    a4 = -158143,
    c1 = 457.8275,
    c2 = 1378661,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14171) "Sr(Ac)2,aq";

  //647
  constant DataRecord Ba_Ac_p1(
    R = Modelica.Constants.R/Ba_Ac_p1.MM,
    G_ref = -935752,
    H_ref = -1016084,
    S_ref = 140.164,
    a1 = 2.772e-05,
    a2 = 3511.46,
    a3 = 0.00010282,
    a4 = -130788,
    c1 = 204.9348,
    c2 = 489507,
    w_ref = 19200,
    z = 1,
    a_0 = 3.72,
    MM = 0.19638) "Ba(Ac)+";

  //648
  constant DataRecord Ba_Ac2_aq(
    R = Modelica.Constants.R/Ba_Ac2_aq.MM,
    G_ref = -1308002,
    H_ref = -1497914,
    S_ref = 249.785,
    a1 = 5.8235e-05,
    a2 = 10964.84,
    a3 = -0.00019061,
    a4 = -161599,
    c1 = 433.6063,
    c2 = 1294479,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19143) "Ba(Ac)2,aq";

  //649
  constant DataRecord Cu_Ac_aq(
    R = Modelica.Constants.R/Cu_Ac_aq.MM,
    G_ref = -321206,
    H_ref = -418274,
    S_ref = 120.081,
    a1 = 3.0547e-05,
    a2 = 4204.21,
    a3 = 7.5086e-05,
    a4 = -133650,
    c1 = 234.3772,
    c2 = 602006,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.090584) "Cu(Ac),aq";

  //650
  constant DataRecord Cu_Ac2_n1(
    R = Modelica.Constants.R/Cu_Ac2_n1.MM,
    G_ref = -690360,
    H_ref = -919392,
    S_ref = 154.808,
    a1 = 6.3059e-05,
    a2 = 12142.18,
    a3 = -0.00023678,
    a4 = -166465,
    c1 = 533.696,
    c2 = 1495106,
    w_ref = 447310,
    z = -1,
    a_0 = 3.72,
    MM = 0.18163) "Cu(Ac)2-";

  //651
  constant DataRecord Rb_Ac_aq(
    R = Modelica.Constants.R/Rb_Ac_aq.MM,
    G_ref = -653164,
    H_ref = -731991,
    S_ref = 224.681,
    a1 = 4.4747e-05,
    a2 = 7671.45,
    a3 = -6.1183e-05,
    a4 = -147984,
    c1 = 142.2397,
    c2 = 281763,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.11251) "Rb(Ac),aq";

  //652
  constant DataRecord Rb_Ac2_n1(
    R = Modelica.Constants.R/Rb_Ac2_n1.MM,
    G_ref = -1020896,
    H_ref = -1223778,
    S_ref = 285.767,
    a1 = 7.8209e-05,
    a2 = 15840.96,
    a3 = -0.00038208,
    a4 = -181757,
    c1 = 335.615,
    c2 = 870272,
    w_ref = 248570,
    z = -1,
    a_0 = 3.72,
    MM = 0.20356) "Rb(Ac)2-";

  //653
  constant DataRecord Tl_Ac_aq(
    R = Modelica.Constants.R/Tl_Ac_aq.MM,
    G_ref = -401078,
    H_ref = -474256,
    S_ref = 230.957,
    a1 = 4.726e-05,
    a2 = 8284.45,
    a3 = -8.5132e-05,
    a4 = -150515,
    c1 = 135.6193,
    c2 = 258751,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.23141) "Tl(Ac),aq";

  //654
  constant DataRecord Tl_Ac2_n1(
    R = Modelica.Constants.R/Tl_Ac2_n1.MM,
    G_ref = -768182,
    H_ref = -964914,
    S_ref = 294.135,
    a1 = 8.0967e-05,
    a2 = 16515.13,
    a3 = -0.00040871,
    a4 = -184544,
    c1 = 321.5483,
    c2 = 825373,
    w_ref = 236100,
    z = -1,
    a_0 = 3.72,
    MM = 0.32246) "Tl(Ac)2-";

  //655
  constant DataRecord Cs_Ac_aq(
    R = Modelica.Constants.R/Cs_Ac_aq.MM,
    G_ref = -661114,
    H_ref = -737723,
    S_ref = 240.58,
    a1 = 4.9315e-05,
    a2 = 8785.31,
    a3 = -0.00010468,
    a4 = -152586,
    c1 = 124.0887,
    c2 = 218673,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15995) "Cs(Ac),aq";

  //656
  constant DataRecord Cs_Ac2_n1(
    R = Modelica.Constants.R/Cs_Ac2_n1.MM,
    G_ref = -1028846,
    H_ref = -1228297,
    S_ref = 306.269,
    a1 = 8.3194e-05,
    a2 = 17056.91,
    a3 = -0.00042958,
    a4 = -186782,
    c1 = 297.3732,
    c2 = 747179,
    w_ref = 217900,
    z = -1,
    a_0 = 3.72,
    MM = 0.251) "Cs(Ac)2-";

  //657
  constant DataRecord Pb_Ac_p1(
    R = Modelica.Constants.R/Pb_Ac_p1.MM,
    G_ref = -406894,
    H_ref = -484842,
    S_ref = 150.624,
    a1 = 2.575e-05,
    a2 = 3030.43,
    a3 = 0.00012168,
    a4 = -128796,
    c1 = 201.1768,
    c2 = 481834,
    w_ref = 2380,
    z = 1,
    a_0 = 3.72,
    MM = 0.26623) "Pb(Ac)+";

  //658
  constant DataRecord Pb_Ac2_aq(
    R = Modelica.Constants.R/Pb_Ac2_aq.MM,
    G_ref = -781948,
    H_ref = -960061,
    S_ref = 292.043,
    a1 = 5.6103e-05,
    a2 = 10442.59,
    a3 = -0.00016975,
    a4 = -159440,
    c1 = 429.3006,
    c2 = 1279513,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.26128) "Pb(Ac)2,aq";

  //659
  constant DataRecord Mn_Ac_p1(
    R = Modelica.Constants.R/Mn_Ac_p1.MM,
    G_ref = -604295,
    H_ref = -709439,
    S_ref = 26.359,
    a1 = 2.5429e-05,
    a2 = 2952.65,
    a3 = 0.00012463,
    a4 = -128474,
    c1 = 265.9635,
    c2 = 646750,
    w_ref = 190580,
    z = 1,
    a_0 = 3.72,
    MM = 0.11398) "Mn(Ac)+";

  //660
  constant DataRecord Mn_Ac2_aq(
    R = Modelica.Constants.R/Mn_Ac2_aq.MM,
    G_ref = -978428,
    H_ref = -1203611,
    S_ref = 101.253,
    a1 = 5.5037e-05,
    a2 = 10184.07,
    a3 = -0.00015998,
    a4 = -158369,
    c1 = 521.8766,
    c2 = 1601284,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.10903) "Mn(Ac)2,aq";

  //661
  constant DataRecord Mn_Ac3_n1(
    R = Modelica.Constants.R/Mn_Ac3_n1.MM,
    G_ref = -1349675,
    H_ref = -1708244,
    S_ref = 131.796,
    a1 = 9.0465e-05,
    a2 = 18833.19,
    a3 = -0.00049961,
    a4 = -194125,
    c1 = 884.0943,
    c2 = 2701680,
    w_ref = 482670,
    z = -1,
    a_0 = 3.72,
    MM = 0.23207) "Mn(Ac)3-";

  //662
  constant DataRecord Co_Ac_p1(
    R = Modelica.Constants.R/Co_Ac_p1.MM,
    G_ref = -432040,
    H_ref = -552623,
    S_ref = -27.196,
    a1 = 2.1043e-05,
    a2 = 1882.47,
    a3 = 0.00016655,
    a4 = -124051,
    c1 = 252.9408,
    c2 = 575798,
    w_ref = 270790,
    z = 1,
    a_0 = 3.72,
    MM = 0.11797) "Co(Ac)+";

  //663
  constant DataRecord Co_Ac2_aq(
    R = Modelica.Constants.R/Co_Ac2_aq.MM,
    G_ref = -806550,
    H_ref = -1052109,
    S_ref = 31.38,
    a1 = 4.9849e-05,
    a2 = 8916.94,
    a3 = -0.00011013,
    a4 = -153130,
    c1 = 482.0474,
    c2 = 1462848,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.11302) "Co(Ac)2,aq";

  //664
  constant DataRecord Co_Ac3_n1(
    R = Modelica.Constants.R/Co_Ac3_n1.MM,
    G_ref = -1179428,
    H_ref = -1563686,
    S_ref = 43.514,
    a1 = 8.5133e-05,
    a2 = 17530.5,
    a3 = -0.00044822,
    a4 = -188740,
    c1 = 829.0261,
    c2 = 2467694,
    w_ref = 615630,
    z = -1,
    a_0 = 3.72,
    MM = 0.23606) "Co(Ac)3-";

  //665
  constant DataRecord Ni_Ac_p1(
    R = Modelica.Constants.R/Ni_Ac_p1.MM,
    G_ref = -423086,
    H_ref = -549987,
    S_ref = -48.953,
    a1 = 1.8224e-05,
    a2 = 1192.94,
    a3 = 0.0001939,
    a4 = -121202,
    c1 = 234.5638,
    c2 = 501013,
    w_ref = 304890,
    z = 1,
    a_0 = 3.72,
    MM = 0.085754) "Ni(Ac)+";

  //666
  constant DataRecord Ni_Ac2_aq(
    R = Modelica.Constants.R/Ni_Ac2_aq.MM,
    G_ref = -797512,
    H_ref = -1051607,
    S_ref = 2.929,
    a1 = 4.6579e-05,
    a2 = 8118.26,
    a3 = -7.8663e-05,
    a4 = -149829,
    c1 = 440.0656,
    c2 = 1316927,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1128) "Ni(Ac)2,aq";

  //667
  constant DataRecord Ni_Ac3_n1(
    R = Modelica.Constants.R/Ni_Ac3_n1.MM,
    G_ref = -1170223,
    H_ref = -1565335,
    S_ref = 7.95,
    a1 = 8.1677e-05,
    a2 = 16686.92,
    a3 = -0.00041516,
    a4 = -185255,
    c1 = 763.1373,
    c2 = 2221064,
    w_ref = 670700,
    z = -1,
    a_0 = 3.72,
    MM = 0.23584) "Ni(Ac)3-";

  //668
  constant DataRecord Cu_Ac_p1(
    R = Modelica.Constants.R/Cu_Ac_p1.MM,
    G_ref = -316478,
    H_ref = -431454,
    S_ref = -5.439,
    a1 = 2.0804e-05,
    a2 = 1825.06,
    a3 = 0.00016857,
    a4 = -123813,
    c1 = 261.4791,
    c2 = 616069,
    w_ref = 237690,
    z = 1,
    a_0 = 3.72,
    MM = 0.12258) "Cu(Ac)+";

  //669
  constant DataRecord Cu_Ac2_aq(
    R = Modelica.Constants.R/Cu_Ac2_aq.MM,
    G_ref = -693791,
    H_ref = -931735,
    S_ref = 59.413,
    a1 = 4.9706e-05,
    a2 = 8881.13,
    a3 = -0.00010847,
    a4 = -152984,
    c1 = 504.6532,
    c2 = 1541419,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.11763) "Cu(Ac)2,aq";

  //670
  constant DataRecord Cu_Ac3_n1(
    R = Modelica.Constants.R/Cu_Ac3_n1.MM,
    G_ref = -1070309,
    H_ref = -1444819,
    S_ref = 79.454,
    a1 = 8.479e-05,
    a2 = 17448.07,
    a3 = -0.00044527,
    a4 = -188401,
    c1 = 862.1969,
    c2 = 2600498,
    w_ref = 560990,
    z = -1,
    a_0 = 3.72,
    MM = 0.24067) "Cu(Ac)3-";

  //671
  constant DataRecord NH4_Ac_aq(
    R = Modelica.Constants.R/NH4_Ac_aq.MM,
    G_ref = -449989,
    H_ref = -616010,
    S_ref = 212.129,
    a1 = 4.7216e-05,
    a2 = 8272.56,
    a3 = -8.4462e-05,
    a4 = -150469,
    c1 = 245.6322,
    c2 = 641127,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.045082) "NH4(Ac),aq";

  //672
  constant DataRecord NH4_Ac2_n1(
    R = Modelica.Constants.R/NH4_Ac2_n1.MM,
    G_ref = -818558,
    H_ref = -1109597,
    S_ref = 270.705,
    a1 = 8.1038e-05,
    a2 = 16530.57,
    a3 = -0.00040893,
    a4 = -184606,
    c1 = 539.4791,
    c2 = 1571431,
    w_ref = 271750,
    z = -1,
    a_0 = 3.72,
    MM = 0.13613) "NH4(Ac)2-";

  //673
  constant DataRecord UO2_Ac_p1(
    R = Modelica.Constants.R/UO2_Ac_p1.MM,
    G_ref = -1339298,
    H_ref = -1521131,
    S_ref = -7.113,
    a1 = 4.0162e-05,
    a2 = 6551.64,
    a3 = -1.715e-05,
    a4 = -143352,
    c1 = 349.4883,
    c2 = 920974,
    w_ref = 240830,
    z = 1,
    a_0 = 3.72,
    MM = 0.32907) "UO2(Ac)+";

  //674
  constant DataRecord UO2_Ac2_aq(
    R = Modelica.Constants.R/UO2_Ac2_aq.MM,
    G_ref = -1723139,
    H_ref = -2027985,
    S_ref = 57.321,
    a1 = 7.1264e-05,
    a2 = 14143.76,
    a3 = -0.00031509,
    a4 = -174741,
    c1 = 675.8114,
    c2 = 2136317,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.29212) "UO2(Ac)2,aq";

  //675
  constant DataRecord UO2_Ac3_n1(
    R = Modelica.Constants.R/UO2_Ac3_n1.MM,
    G_ref = -2102042,
    H_ref = -2543705,
    S_ref = 76.567,
    a1 = 0.00010883,
    a2 = 23317.18,
    a3 = -0.00067601,
    a4 = -212664,
    c1 = 1151.9435,
    c2 = 3606001,
    w_ref = 565930,
    z = -1,
    a_0 = 3.72,
    MM = 0.44716) "UO2(Ac)3-";

  //676
  constant DataRecord Ag_Ac_aq(
    R = Modelica.Constants.R/Ag_Ac_aq.MM,
    G_ref = -296395,
    H_ref = -383464,
    S_ref = 162.758,
    a1 = 3.514e-05,
    a2 = 5324.43,
    a3 = 3.1317e-05,
    a4 = -138281,
    c1 = 202.3776,
    c2 = 490783,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13491) "Ag(Ac),aq";

  //677
  constant DataRecord Ag_Ac2_n1(
    R = Modelica.Constants.R/Ag_Ac2_n1.MM,
    G_ref = -665214,
    H_ref = -880983,
    S_ref = 208.782,
    a1 = 6.7901e-05,
    a2 = 13322.61,
    a3 = -0.00028281,
    a4 = -171343,
    c1 = 463.7433,
    c2 = 1278099,
    w_ref = 365720,
    z = -1,
    a_0 = 3.72,
    MM = 0.22596) "Ag(Ac)2-";

  //678
  constant DataRecord Cd_Ac_p1(
    R = Modelica.Constants.R/Cd_Ac_p1.MM,
    G_ref = -457981,
    H_ref = -568689,
    S_ref = 27.614,
    a1 = 2.6377e-05,
    a2 = 3184.53,
    a3 = 0.00011545,
    a4 = -129432,
    c1 = 269.0458,
    c2 = 658256,
    w_ref = 188110,
    z = 1,
    a_0 = 3.72,
    MM = 0.17144) "Cd(Ac)+";

  //679
  constant DataRecord Cd_Ac2_aq(
    R = Modelica.Constants.R/Cd_Ac2_aq.MM,
    G_ref = -834290,
    H_ref = -1064912,
    S_ref = 102.926,
    a1 = 5.6103e-05,
    a2 = 10442.59,
    a3 = -0.00016975,
    a4 = -159440,
    c1 = 528.3354,
    c2 = 1623731,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16649) "Cd(Ac)2,aq";

  //680
  constant DataRecord Cd_Ac3_n1(
    R = Modelica.Constants.R/Cd_Ac3_n1.MM,
    G_ref = -1205118,
    H_ref = -1569335,
    S_ref = 133.47,
    a1 = 9.1643e-05,
    a2 = 19122.39,
    a3 = -0.00051128,
    a4 = -195322,
    c1 = 894.7505,
    c2 = 2739625,
    w_ref = 479860,
    z = -1,
    a_0 = 3.72,
    MM = 0.28953) "Cd(Ac)3-";

  //681
  constant DataRecord Hg_Ac_p1(
    R = Modelica.Constants.R/Hg_Ac_p1.MM,
    G_ref = -229116,
    H_ref = -332168,
    S_ref = 77.404,
    a1 = 2.3573e-05,
    a2 = 2501.49,
    a3 = 0.00014198,
    a4 = -126612,
    c1 = 293.6181,
    c2 = 767563,
    w_ref = 113470,
    z = 1,
    a_0 = 3.72,
    MM = 0.25963) "Hg(Ac)+";

  //682
  constant DataRecord Hg_Ac2_aq(
    R = Modelica.Constants.R/Hg_Ac2_aq.MM,
    G_ref = -613291,
    H_ref = -831696,
    S_ref = 167.778,
    a1 = 5.326e-05,
    a2 = 9748.89,
    a3 = -0.00014259,
    a4 = -156569,
    c1 = 589.6938,
    c2 = 1836998,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.25468) "Hg(Ac)2,aq";

  //683
  constant DataRecord Hg_Ac3_n1(
    R = Modelica.Constants.R/Hg_Ac3_n1.MM,
    G_ref = -1000060,
    H_ref = -1346830,
    S_ref = 215.058,
    a1 = 8.8057e-05,
    a2 = 18244.79,
    a3 = -0.00047639,
    a4 = -191694,
    c1 = 987.0458,
    c2 = 3100089,
    w_ref = 355970,
    z = -1,
    a_0 = 3.72,
    MM = 0.37772) "Hg(Ac)3-";

  //684
  constant DataRecord Sc_Ac_p2(
    R = Modelica.Constants.R/Sc_Ac_p2.MM,
    G_ref = -974914,
    H_ref = -1121730,
    S_ref = -177.402,
    a1 = 1.137e-05,
    a2 = -478.52,
    a3 = 0.00025914,
    a4 = -114290,
    c1 = 252.3768,
    c2 = 432617,
    w_ref = 711820,
    z = 2,
    a_0 = 3.72,
    MM = 0.104) "Sc(Ac)+2";

  //685
  constant DataRecord Sc_Ac2_p1(
    R = Modelica.Constants.R/Sc_Ac2_p1.MM,
    G_ref = -1358294,
    H_ref = -1628915,
    S_ref = -115.06,
    a1 = 3.8825e-05,
    a2 = 6223.16,
    a3 = -3.761e-06,
    a4 = -141997,
    c1 = 445.6726,
    c2 = 1202356,
    w_ref = 406100,
    z = 1,
    a_0 = 3.72,
    MM = 0.16305) "Sc(Ac)2+";

  //686
  constant DataRecord Sc_Ac3_aq(
    R = Modelica.Constants.R/Sc_Ac3_aq.MM,
    G_ref = -1737908,
    H_ref = -2141539,
    S_ref = -83.68,
    a1 = 6.9152e-05,
    a2 = 13629.3,
    a3 = -0.00029514,
    a4 = -172611,
    c1 = 602.2224,
    c2 = 1880541,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.12609) "Sc(Ac)3,aq";

  //687
  constant DataRecord U_Ac_p2(
    R = Modelica.Constants.R/U_Ac_p2.MM,
    G_ref = -859728,
    H_ref = -985165,
    S_ref = -89.119,
    a1 = 2.036e-05,
    a2 = 1716.74,
    a3 = 0.00017283,
    a4 = -123365,
    c1 = 433.7381,
    c2 = 1105722,
    w_ref = 578350,
    z = 2,
    a_0 = 3.72,
    MM = 0.29707) "U(Ac)+2";

  //688
  constant DataRecord U_Ac2_p1(
    R = Modelica.Constants.R/U_Ac2_p1.MM,
    G_ref = -1242439,
    H_ref = -1482391,
    S_ref = 5.439,
    a1 = 4.8725e-05,
    a2 = 8642.22,
    a3 = -9.9236e-05,
    a4 = -151996,
    c1 = 806.6317,
    c2 = 2515663,
    w_ref = 222760,
    z = 1,
    a_0 = 3.72,
    MM = 0.35612) "U(Ac)2+";

  //689
  constant DataRecord U_Ac3_aq(
    R = Modelica.Constants.R/U_Ac3_aq.MM,
    G_ref = -1620714,
    H_ref = -1982337,
    S_ref = 74.057,
    a1 = 8.0873e-05,
    a2 = 16490.32,
    a3 = -0.00040741,
    a4 = -184439,
    c1 = 1240.8652,
    c2 = 4100303,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.31916) "U(Ac)3,aq";

  //690
  constant DataRecord Pr_Ac_p2(
    R = Modelica.Constants.R/Pr_Ac_p2.MM,
    G_ref = -1065121,
    H_ref = -1204490,
    S_ref = -112.131,
    a1 = 1.0911e-05,
    a2 = -590.61,
    a3 = 0.00026358,
    a4 = -113826,
    c1 = 175.487,
    c2 = 196748,
    w_ref = 613830,
    z = 2,
    a_0 = 3.72,
    MM = 0.19995) "Pr(Ac)+2";

  //691
  constant DataRecord Pr_Ac2_p1(
    R = Modelica.Constants.R/Pr_Ac2_p1.MM,
    G_ref = -1445572,
    H_ref = -1701675,
    S_ref = -25.941,
    a1 = 3.8226e-05,
    a2 = 6077.26,
    a3 = 1.845e-06,
    a4 = -141394,
    c1 = 300.8007,
    c2 = 742150,
    w_ref = 270790,
    z = 1,
    a_0 = 3.72,
    MM = 0.259) "Pr(Ac)2+";

  //692
  constant DataRecord Pr_Ac3_aq(
    R = Modelica.Constants.R/Pr_Ac3_aq.MM,
    G_ref = -1822927,
    H_ref = -2203922,
    S_ref = 32.635,
    a1 = 6.8993e-05,
    a2 = 13589.88,
    a3 = -0.00029341,
    a4 = -172448,
    c1 = 378.4302,
    c2 = 1102702,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22204) "Pr(Ac)3,aq";

  //693
  constant DataRecord La_Ac_p2(
    R = Modelica.Constants.R/La_Ac_p2.MM,
    G_ref = -1070058,
    H_ref = -1207963,
    S_ref = -123.846,
    a1 = 1.32e-05,
    a2 = -33.64,
    a3 = 0.0002421,
    a4 = -116131,
    c1 = 234.9442,
    c2 = 398099,
    w_ref = 630400,
    z = 2,
    a_0 = 3.72,
    MM = 0.19795) "La(Ac)+2";

  //694
  constant DataRecord La_Ac2_p1(
    R = Modelica.Constants.R/La_Ac2_p1.MM,
    G_ref = -1448333,
    H_ref = -1704269,
    S_ref = -42.258,
    a1 = 4.0789e-05,
    a2 = 6703.44,
    a3 = -2.2861e-05,
    a4 = -143980,
    c1 = 415.8766,
    c2 = 1135006,
    w_ref = 293010,
    z = 1,
    a_0 = 3.72,
    MM = 0.257) "La(Ac)2+";

  //695
  constant DataRecord La_Ac3_aq(
    R = Modelica.Constants.R/La_Ac3_aq.MM,
    G_ref = -1826525,
    H_ref = -2209093,
    S_ref = 11.715,
    a1 = 7.1765e-05,
    a2 = 14268.28,
    a3 = -0.00032042,
    a4 = -175255,
    c1 = 569.4725,
    c2 = 1766711,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22004) "La(Ac)3,aq";

  //696
  constant DataRecord Nd_Ac_p2(
    R = Modelica.Constants.R/Nd_Ac_p2.MM,
    G_ref = -1056502,
    H_ref = -1194406,
    S_ref = -109.202,
    a1 = 1.026e-05,
    a2 = -751.49,
    a3 = 0.00027026,
    a4 = -113160,
    c1 = 199.9404,
    c2 = 283039,
    w_ref = 609780,
    z = 2,
    a_0 = 3.72,
    MM = 0.20328) "Nd(Ac)+2";

  //697
  constant DataRecord Nd_Ac2_p1(
    R = Modelica.Constants.R/Nd_Ac2_p1.MM,
    G_ref = -1436493,
    H_ref = -1690796,
    S_ref = -22.175,
    a1 = 3.7491e-05,
    a2 = 5899.02,
    a3 = 8.644e-06,
    a4 = -140658,
    c1 = 348.5975,
    c2 = 910518,
    w_ref = 263800,
    z = 1,
    a_0 = 3.72,
    MM = 0.26233) "Nd(Ac)2+";

  //698
  constant DataRecord Nd_Ac3_aq(
    R = Modelica.Constants.R/Nd_Ac3_aq.MM,
    G_ref = -1814015,
    H_ref = -2192793,
    S_ref = 38.074,
    a1 = 6.8202e-05,
    a2 = 13397.84,
    a3 = -0.00028616,
    a4 = -171657,
    c1 = 460.3053,
    c2 = 1387276,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22537) "Nd(Ac)3,aq";

  //699
  constant DataRecord Ce_Ac_p2(
    R = Modelica.Constants.R/Ce_Ac_p2.MM,
    G_ref = -1061021,
    H_ref = -1198256,
    S_ref = -106.274,
    a1 = 1.2337e-05,
    a2 = -243.72,
    a3 = 0.00025024,
    a4 = -115261,
    c1 = 224.5791,
    c2 = 371251,
    w_ref = 601740,
    z = 2,
    a_0 = 3.72,
    MM = 0.19916) "Ce(Ac)+2";

  //700
  constant DataRecord Ce_Ac2_p1(
    R = Modelica.Constants.R/Ce_Ac2_p1.MM,
    G_ref = -1441597,
    H_ref = -1694897,
    S_ref = -17.991,
    a1 = 3.9814e-05,
    a2 = 6465.41,
    a3 = -1.3468e-05,
    a4 = -142997,
    c1 = 397.4892,
    c2 = 1082627,
    w_ref = 257020,
    z = 1,
    a_0 = 3.72,
    MM = 0.25821) "Ce(Ac)2+";

  //701
  constant DataRecord Ce_Ac3_aq(
    R = Modelica.Constants.R/Ce_Ac3_aq.MM,
    G_ref = -1819036,
    H_ref = -2196433,
    S_ref = 43.095,
    a1 = 7.0815e-05,
    a2 = 14034.22,
    a3 = -0.0003108,
    a4 = -174289,
    c1 = 543.9999,
    c2 = 1678177,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22125) "Ce(Ac)3,aq";

  //702
  constant DataRecord Gd_Ac_p2(
    R = Modelica.Constants.R/Gd_Ac_p2.MM,
    G_ref = -1048050,
    H_ref = -1184490,
    S_ref = -107.529,
    a1 = 1.1968e-05,
    a2 = -332.42,
    a3 = 0.00025341,
    a4 = -114897,
    c1 = 239.844,
    c2 = 423028,
    w_ref = 605720,
    z = 2,
    a_0 = 3.72,
    MM = 0.21629) "Gd(Ac)+2";

  //703
  constant DataRecord Gd_Ac2_p1(
    R = Modelica.Constants.R/Gd_Ac2_p1.MM,
    G_ref = -1428208,
    H_ref = -1680880,
    S_ref = -19.665,
    a1 = 3.9399e-05,
    a2 = 6365.29,
    a3 = -9.799e-06,
    a4 = -142582,
    c1 = 426.8642,
    c2 = 1183649,
    w_ref = 260370,
    z = 1,
    a_0 = 3.72,
    MM = 0.27534) "Gd(Ac)2+";

  //704
  constant DataRecord Gd_Ac3_aq(
    R = Modelica.Constants.R/Gd_Ac3_aq.MM,
    G_ref = -1805354,
    H_ref = -2182291,
    S_ref = 41.003,
    a1 = 7.034e-05,
    a2 = 13918.58,
    a3 = -0.00030633,
    a4 = -173808,
    c1 = 593.1251,
    c2 = 1848922,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.23838) "Gd(Ac)3,aq";

  //705
  constant DataRecord Sm_Ac_p2(
    R = Modelica.Constants.R/Sm_Ac_p2.MM,
    G_ref = -1051188,
    H_ref = -1190557,
    S_ref = -116.315,
    a1 = 1.0989e-05,
    a2 = -571.83,
    a3 = 0.00026287,
    a4 = -113905,
    c1 = 200.14,
    c2 = 281123,
    w_ref = 617930,
    z = 2,
    a_0 = 3.72,
    MM = 0.20939) "Sm(Ac)+2";

  //706
  constant DataRecord Sm_Ac2_p1(
    R = Modelica.Constants.R/Sm_Ac2_p1.MM,
    G_ref = -1431723,
    H_ref = -1688244,
    S_ref = -31.798,
    a1 = 3.8321e-05,
    a2 = 6101.9,
    a3 = 5.77e-07,
    a4 = -141495,
    c1 = 348.8272,
    c2 = 906777,
    w_ref = 277980,
    z = 1,
    a_0 = 3.72,
    MM = 0.26844) "Sm(Ac)2+";

  //707
  constant DataRecord Sm_Ac3_aq(
    R = Modelica.Constants.R/Sm_Ac3_aq.MM,
    G_ref = -1810124,
    H_ref = -2192454,
    S_ref = 25.104,
    a1 = 6.9073e-05,
    a2 = 13610.84,
    a3 = -0.0002946,
    a4 = -172536,
    c1 = 458.4856,
    c2 = 1380954,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.23148) "Sm(Ac)3,aq";

  //708
  constant DataRecord Yb_Ac_p2(
    R = Modelica.Constants.R/Yb_Ac_p2.MM,
    G_ref = -1024076,
    H_ref = -1171687,
    S_ref = -153.134,
    a1 = 9.5839e-06,
    a2 = -916.51,
    a3 = 0.00027679,
    a4 = -112478,
    c1 = 243.3908,
    c2 = 413438,
    w_ref = 674170,
    z = 2,
    a_0 = 3.72,
    MM = 0.23208) "Yb(Ac)+2";

  //709
  constant DataRecord Yb_Ac2_p1(
    R = Modelica.Constants.R/Yb_Ac2_p1.MM,
    G_ref = -1403816,
    H_ref = -1672554,
    S_ref = -82.006,
    a1 = 3.6799e-05,
    a2 = 5730.36,
    a3 = 1.518e-05,
    a4 = -139959,
    c1 = 430.06,
    c2 = 1164939,
    w_ref = 353510,
    z = 1,
    a_0 = 3.72,
    MM = 0.29113) "Yb(Ac)2+";

  //710
  constant DataRecord Yb_Ac3_aq(
    R = Modelica.Constants.R/Yb_Ac3_aq.MM,
    G_ref = -1780669,
    H_ref = -2179404,
    S_ref = -40.585,
    a1 = 6.7093e-05,
    a2 = 13127.01,
    a3 = -0.00027549,
    a4 = -170536,
    c1 = 584.0274,
    c2 = 1817304,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.25417) "Yb(Ac)3,aq";

  //711
  constant DataRecord Eu_Ac_p2(
    R = Modelica.Constants.R/Eu_Ac_p2.MM,
    G_ref = -959768,
    H_ref = -1105748,
    S_ref = -130.122,
    a1 = 1.1506e-05,
    a2 = -446.27,
    a3 = 0.00025811,
    a4 = -114424,
    c1 = 239.0357,
    c2 = 409605,
    w_ref = 638850,
    z = 2,
    a_0 = 3.72,
    MM = 0.211) "Eu(Ac)+2";

  //712
  constant DataRecord Eu_Ac2_p1(
    R = Modelica.Constants.R/Eu_Ac2_p1.MM,
    G_ref = -1340637,
    H_ref = -1605275,
    S_ref = -50.208,
    a1 = 3.8923e-05,
    a2 = 6247,
    a3 = -4.699e-06,
    a4 = -142093,
    c1 = 423.8045,
    c2 = 1157458,
    w_ref = 308950,
    z = 1,
    a_0 = 3.72,
    MM = 0.27005) "Eu(Ac)2+";

  //713
  constant DataRecord Eu_Ac3_aq(
    R = Modelica.Constants.R/Eu_Ac3_aq.MM,
    G_ref = -1718327,
    H_ref = -2110075,
    S_ref = 0.837,
    a1 = 6.9627e-05,
    a2 = 13744.94,
    a3 = -0.00029959,
    a4 = -173092,
    c1 = 580.3886,
    c2 = 1804655,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.23309) "Eu(Ac)3,aq";

  //714
  constant DataRecord Y_Ac_p2(
    R = Modelica.Constants.R/Y_Ac_p2.MM,
    G_ref = -1069723,
    H_ref = -1218088,
    S_ref = -171.544,
    a1 = 1.2039e-05,
    a2 = -316.98,
    a3 = 0.0002532,
    a4 = -114960,
    c1 = 249.8325,
    c2 = 426864,
    w_ref = 702160,
    z = 2,
    a_0 = 3.72,
    MM = 0.14795) "Y(Ac)+2";

  //715
  constant DataRecord Y_Ac2_p1(
    R = Modelica.Constants.R/Y_Ac2_p1.MM,
    G_ref = -1449840,
    H_ref = -1721381,
    S_ref = -107.11,
    a1 = 3.9569e-05,
    a2 = 6406.83,
    a3 = -1.1447e-05,
    a4 = -142754,
    c1 = 441.4045,
    c2 = 1191130,
    w_ref = 394840,
    z = 1,
    a_0 = 3.72,
    MM = 0.207) "Y(Ac)2+";

  //716
  constant DataRecord Y_Ac3_aq(
    R = Modelica.Constants.R/Y_Ac3_aq.MM,
    G_ref = -1826944,
    H_ref = -2230783,
    S_ref = -73.22,
    a1 = 7.0023e-05,
    a2 = 13842.26,
    a3 = -0.00030355,
    a4 = -173494,
    c1 = 596.7639,
    c2 = 1861570,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17004) "Y(Ac)3,aq";

  //717
  constant DataRecord Lu_Ac_p2(
    R = Modelica.Constants.R/Lu_Ac_p2.MM,
    G_ref = -1051397,
    H_ref = -1207226,
    S_ref = -189.954,
    a1 = 9.4592e-06,
    a2 = -946.5,
    a3 = 0.00027787,
    a4 = -112357,
    c1 = 272.9583,
    c2 = 497817,
    w_ref = 731610,
    z = 2,
    a_0 = 3.72,
    MM = 0.234) "Lu(Ac)+2";

  //718
  constant DataRecord Lu_Ac2_p1(
    R = Modelica.Constants.R/Lu_Ac2_p1.MM,
    G_ref = -1431556,
    H_ref = -1712553,
    S_ref = -132.214,
    a1 = 3.6702e-05,
    a2 = 5706.22,
    a3 = 1.623e-05,
    a4 = -139859,
    c1 = 484.4679,
    c2 = 1329566,
    w_ref = 429950,
    z = 1,
    a_0 = 3.72,
    MM = 0.29305) "Lu(Ac)2+";

  //719
  constant DataRecord Lu_Ac3_aq(
    R = Modelica.Constants.R/Lu_Ac3_aq.MM,
    G_ref = -1808701,
    H_ref = -2224298,
    S_ref = -105.855,
    a1 = 6.6697e-05,
    a2 = 13029.77,
    a3 = -0.00027153,
    a4 = -170134,
    c1 = 664.0832,
    c2 = 2095556,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.25609) "Lu(Ac)3,aq";

  //720
  constant DataRecord Tb_Ac_p2(
    R = Modelica.Constants.R/Tb_Ac_p2.MM,
    G_ref = -1051816,
    H_ref = -1198298,
    S_ref = -135.98,
    a1 = 1.2237e-05,
    a2 = -267.02,
    a3 = 0.00025089,
    a4 = -115165,
    c1 = 218.3123,
    c2 = 334816,
    w_ref = 647470,
    z = 2,
    a_0 = 3.72,
    MM = 0.21796) "Tb(Ac)+2";

  //721
  constant DataRecord Tb_Ac2_p1(
    R = Modelica.Constants.R/Tb_Ac2_p1.MM,
    G_ref = -1431974,
    H_ref = -1697784,
    S_ref = -58.576,
    a1 = 3.9733e-05,
    a2 = 6447.29,
    a3 = -1.3104e-05,
    a4 = -142921,
    c1 = 382.5921,
    c2 = 1011536,
    w_ref = 317310,
    z = 1,
    a_0 = 3.72,
    MM = 0.27701) "Tb(Ac)2+";

  //722
  constant DataRecord Tb_Ac3_aq(
    R = Modelica.Constants.R/Tb_Ac3_aq.MM,
    G_ref = -1809120,
    H_ref = -2202750,
    S_ref = -9.623,
    a1 = 7.0498e-05,
    a2 = 13957.95,
    a3 = -0.00030804,
    a4 = -173971,
    c1 = 509.4309,
    c2 = 1558021,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.24005) "Tb(Ac)3,aq";

  //723
  constant DataRecord Pb_Ac3_n1(
    R = Modelica.Constants.R/Pb_Ac3_n1.MM,
    G_ref = -1162106,
    H_ref = -1459212,
    S_ref = 370.702,
    a1 = 9.0428e-05,
    a2 = 18824.69,
    a3 = -0.00049941,
    a4 = -194092,
    c1 = 694.2227,
    c2 = 2157823,
    w_ref = 120120,
    z = -1,
    a_0 = 3.72,
    MM = 0.38432) "Pb(Ac)3-";

  //724
  constant DataRecord Fe_Ac_p1(
    R = Modelica.Constants.R/Fe_Ac_p1.MM,
    G_ref = -468190,
    H_ref = -581827,
    S_ref = -6.276,
    a1 = 2.186e-05,
    a2 = 2083,
    a3 = 0.00015842,
    a4 = -124884,
    c1 = 247.3413,
    c2 = 565940,
    w_ref = 240830,
    z = 1,
    a_0 = 3.72,
    MM = 0.11489) "Fe(Ac)+";

  //725
  constant DataRecord Fe_Ac2_aq(
    R = Modelica.Constants.R/Fe_Ac2_aq.MM,
    G_ref = -844331,
    H_ref = -1084074,
    S_ref = 48.242,
    a1 = 5.0918e-05,
    a2 = 9178.44,
    a3 = -0.00012046,
    a4 = -154214,
    c1 = 476.0099,
    c2 = 1442936,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.10994) "Fe(Ac)2,aq";

  //726
  constant DataRecord Zn_Ac_p1(
    R = Modelica.Constants.R/Zn_Ac_p1.MM,
    G_ref = -525761,
    H_ref = -649022,
    S_ref = 39.33,
    a1 = 2.0286e-05,
    a2 = 1698.7,
    a3 = 0.00017352,
    a4 = -123294,
    c1 = 252.9755,
    c2 = 607701,
    w_ref = 171540,
    z = 1,
    a_0 = 3.72,
    MM = 0.12441) "Zn(Ac)+";

  //727
  constant DataRecord Zn_Ac2_aq(
    R = Modelica.Constants.R/Zn_Ac2_aq.MM,
    G_ref = -905627,
    H_ref = -1135956,
    S_ref = 94.014,
    a1 = 4.9138e-05,
    a2 = 8743.64,
    a3 = -0.00010337,
    a4 = -152419,
    c1 = 498.3236,
    c2 = 1520495,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.11946) "Zn(Ac)2,aq";

  //728
  constant DataRecord Zn_Ac3_n1(
    R = Modelica.Constants.R/Zn_Ac3_n1.MM,
    G_ref = -1279216,
    H_ref = -1650764,
    S_ref = 104.6,
    a1 = 8.3819e-05,
    a2 = 17211.85,
    a3 = -0.00043621,
    a4 = -187426,
    c1 = 850.1164,
    c2 = 2570503,
    w_ref = 523540,
    z = -1,
    a_0 = 3.72,
    MM = 0.2425) "Zn(Ac)3-";

  //729
  constant DataRecord Ca_Ac_p1(
    R = Modelica.Constants.R/Ca_Ac_p1.MM,
    G_ref = -927425,
    H_ref = -1027674,
    S_ref = 52.3,
    a1 = 2.4686e-05,
    a2 = 2771.15,
    a3 = 0.00013182,
    a4 = -127725,
    c1 = 243.4988,
    c2 = 580978,
    w_ref = 152130,
    z = 1,
    a_0 = 3.72,
    MM = 0.099124) "Ca(Ac)+";

  //730
  constant DataRecord Ca_Ac2_aq(
    R = Modelica.Constants.R/Ca_Ac2_aq.MM,
    G_ref = -1299425,
    H_ref = -1521804,
    S_ref = 135.143,
    a1 = 5.4355e-05,
    a2 = 10015.62,
    a3 = -0.00015295,
    a4 = -157674,
    c1 = 484.9541,
    c2 = 1472948,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.094167) "Ca(Ac)2,aq";

  //731
  constant DataRecord Al_Ac_p2(
    R = Modelica.Constants.R/Al_Ac_p2.MM,
    G_ref = -872029,
    H_ref = -1043322,
    S_ref = -264.01,
    a1 = 1.0214e-05,
    a2 = -762.95,
    a3 = 0.0002708,
    a4 = -113114,
    c1 = 280.333,
    c2 = 488227,
    w_ref = 841610,
    z = 2,
    a_0 = 3.72,
    MM = 0.086024) "Al(Ac)+2";

  //732
  constant DataRecord Al_Ac2_p1(
    R = Modelica.Constants.R/Al_Ac2_p1.MM,
    G_ref = -1256204,
    H_ref = -1560674,
    S_ref = -233.049,
    a1 = 3.7644e-05,
    a2 = 5934.75,
    a3 = 7.552e-06,
    a4 = -140804,
    c1 = 493.1262,
    c2 = 1310860,
    w_ref = 582330,
    z = 1,
    a_0 = 3.72,
    MM = 0.14507) "Al(Ac)2+";

  //733
  constant DataRecord Al_Ac3_aq(
    R = Modelica.Constants.R/Al_Ac3_aq.MM,
    G_ref = -1636362,
    H_ref = -2084804,
    S_ref = -238.488,
    a1 = 6.7172e-05,
    a2 = 13145.46,
    a3 = -0.00027602,
    a4 = -170611,
    c1 = 654.9855,
    c2 = 2063938,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.10811) "Al(Ac)3,aq";

  //734
  constant DataRecord ACETALDEHYDE_aq(
    R = Modelica.Constants.R/ACETALDEHYDE_aq.MM,
    G_ref = -140038,
    H_ref = -210790,
    S_ref = 138.072,
    a1 = 3.1883e-05,
    a2 = 4413.24,
    a3 = 9.2027e-05,
    a4 = -134516,
    c1 = 144.6752,
    c2 = -63296,
    w_ref = -154220,
    z = 0,
    a_0 = 3.72,
    MM = 0.028052) "ACETALDEHYDE,aq";

  //735
  constant DataRecord BUTANAL_aq(
    R = Modelica.Constants.R/BUTANAL_aq.MM,
    G_ref = -120164,
    H_ref = -255517,
    S_ref = 194.974,
    a1 = 5.0101e-05,
    a2 = 8541.72,
    a3 = -1.644e-06,
    a4 = -151582,
    c1 = 289.8851,
    c2 = 111972,
    w_ref = -117070,
    z = 0,
    a_0 = 3.72,
    MM = 0.056103) "BUTANAL,aq";

  //736
  constant DataRecord DECANAL_aq(
    R = Modelica.Constants.R/DECANAL_aq.MM,
    G_ref = -69454,
    H_ref = -398693,
    S_ref = 362.334,
    a1 = 0.00010475,
    a2 = 20929.79,
    a3 = -0.00028329,
    a4 = 202794,
    c1 = 725.5115,
    c2 = 637792,
    w_ref = -5560,
    z = 0,
    a_0 = 3.72,
    MM = 0.14026) "DECANAL,aq";

  //737
  constant DataRecord FORMALDEHYDE_aq(
    R = Modelica.Constants.R/FORMALDEHYDE_aq.MM,
    G_ref = -109328,
    H_ref = -145030,
    S_ref = 119.244,
    a1 = 2.2222e-05,
    a2 = 2223.34,
    a3 = 0.00014184,
    a4 = -125461,
    c1 = 79.3207,
    c2 = -142662,
    w_ref = -166690,
    z = 0,
    a_0 = 3.72,
    MM = 0.014026) "FORMALDEHYDE,aq";

  //738
  constant DataRecord HEPTANAL_aq(
    R = Modelica.Constants.R/HEPTANAL_aq.MM,
    G_ref = -89914,
    H_ref = -322210,
    S_ref = 278.236,
    a1 = 7.7428e-05,
    a2 = 14735.71,
    a3 = -0.00014247,
    a4 = -177188,
    c1 = 507.6983,
    c2 = 374882,
    w_ref = -61300,
    z = 0,
    a_0 = 3.72,
    MM = 0.098181) "HEPTANAL,aq";

  //739
  constant DataRecord HEXANAL_aq(
    R = Modelica.Constants.R/HEXANAL_aq.MM,
    G_ref = -103972,
    H_ref = -303968,
    S_ref = 250.203,
    a1 = 6.8319e-05,
    a2 = 12670.16,
    a3 = -9.5303e-05,
    a4 = -168649,
    c1 = 435.0933,
    c2 = 287248,
    w_ref = -79870,
    z = 0,
    a_0 = 3.72,
    MM = 0.084155) "HEXANAL,aq";

  //740
  constant DataRecord NONANAL_aq(
    R = Modelica.Constants.R/NONANAL_aq.MM,
    G_ref = -75730,
    H_ref = -372627,
    S_ref = 334.302,
    a1 = 9.5646e-05,
    a2 = 18864.28,
    a3 = -0.00023614,
    a4 = -194255,
    c1 = 652.9073,
    c2 = 550154,
    w_ref = -24140,
    z = 0,
    a_0 = 3.72,
    MM = 0.12623) "NONANAL,aq";

  //741
  constant DataRecord OCTANAL_aq(
    R = Modelica.Constants.R/OCTANAL_aq.MM,
    G_ref = -84977,
    H_ref = -349573,
    S_ref = 306.269,
    a1 = 8.6537e-05,
    a2 = 16798.72,
    a3 = -0.00018897,
    a4 = -185715,
    c1 = 580.3024,
    c2 = 922760,
    w_ref = -42720,
    z = 0,
    a_0 = 3.72,
    MM = 0.11221) "OCTANAL,aq";

  //742
  constant DataRecord PENTANAL_aq(
    R = Modelica.Constants.R/PENTANAL_aq.MM,
    G_ref = -113052,
    H_ref = -280746,
    S_ref = 222.17,
    a1 = 5.921e-05,
    a2 = 10604.64,
    a3 = -4.8149e-05,
    a4 = -160109,
    c1 = 362.4892,
    c2 = 199610,
    w_ref = -98490,
    z = 0,
    a_0 = 3.72,
    MM = 0.070129) "PENTANAL,aq";

  //743
  constant DataRecord PROPANAL_aq(
    R = Modelica.Constants.R/PROPANAL_aq.MM,
    G_ref = -136942,
    H_ref = -239994,
    S_ref = 166.105,
    a1 = 4.0992e-05,
    a2 = 6476.16,
    a3 = 4.5522e-05,
    a4 = -143043,
    c1 = 217.2801,
    c2 = 24338,
    w_ref = -135650,
    z = 0,
    a_0 = 3.72,
    MM = 0.042077) "PROPANAL,aq";

  //744
  constant DataRecord UOH_p2(
    R = Modelica.Constants.R/UOH_p2.MM,
    G_ref = -676971,
    H_ref = -701657,
    S_ref = 5.021,
    a1 = 1.1772e-05,
    a2 = -380.2,
    a3 = 0.00025524,
    a4 = -114696,
    c1 = -28.8541,
    c2 = -455943,
    w_ref = 434130,
    z = 2,
    a_0 = 3.72,
    MM = 0.25504) "UOH+2";

  //745
  constant DataRecord UO_p1(
    R = Modelica.Constants.R/UO_p1.MM,
    G_ref = -639734,
    H_ref = -643918,
    S_ref = 73.22,
    a1 = 1.22e-05,
    a2 = -278.15,
    a3 = 0.00025174,
    a4 = -115119,
    c1 = -159.2912,
    c2 = -808788,
    w_ref = 120210,
    z = 1,
    a_0 = 3.72,
    MM = 0.25403) "UO+";

  //746
  constant DataRecord HUO2_aq(
    R = Modelica.Constants.R/HUO2_aq.MM,
    G_ref = -828432,
    H_ref = -858138,
    S_ref = 221.334,
    a1 = 2.1599e-05,
    a2 = 2018.86,
    a3 = 0.00016104,
    a4 = -124616,
    c1 = -329.9264,
    c2 = -1359361,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.27003) "HUO2,aq";

  //747
  constant DataRecord U_OH_p3(
    R = Modelica.Constants.R/U_OH_p3.MM,
    G_ref = -763998,
    H_ref = -830106,
    S_ref = -199.995,
    a1 = 9.042e-06,
    a2 = -1048.68,
    a3 = 0.00028196,
    a4 = -111935,
    c1 = 157.1025,
    c2 = 20476,
    w_ref = 964750,
    z = 3,
    a_0 = 3.72,
    MM = 0.25504) "U(OH)+3";

  //748
  constant DataRecord UO_p2(
    R = Modelica.Constants.R/UO_p2.MM,
    G_ref = -755630,
    H_ref = -830106,
    S_ref = -139.746,
    a1 = 3.3045e-06,
    a2 = -2446.09,
    a3 = 0.00033613,
    a4 = -106156,
    c1 = 30.5963,
    c2 = -320432,
    w_ref = 656220,
    z = 2,
    a_0 = 3.72,
    MM = 0.25403) "UO+2";

  //749
  constant DataRecord HUO2_p1(
    R = Modelica.Constants.R/HUO2_p1.MM,
    G_ref = -975709,
    H_ref = -1066083,
    S_ref = -47.698,
    a1 = 1.3511e-05,
    a2 = 44.14,
    a3 = 0.00023865,
    a4 = -116453,
    c1 = 58.7923,
    c2 = -109918,
    w_ref = 304890,
    z = 1,
    a_0 = 3.72,
    MM = 0.27104) "HUO2+";

  //750
  constant DataRecord UO2_aq(
    R = Modelica.Constants.R/UO2_aq.MM,
    G_ref = -978219,
    H_ref = -1086585,
    S_ref = -108.784,
    a1 = 1.0319e-05,
    a2 = -735.84,
    a3 = 0.00026947,
    a4 = -113227,
    c1 = 71.4795,
    c2 = 35819,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.23803) "UO2,aq";

  //751
  constant DataRecord HUO3_n1(
    R = Modelica.Constants.R/HUO3_n1.MM,
    G_ref = -1146834,
    H_ref = -1343901,
    S_ref = -172.799,
    a1 = 1.3315e-05,
    a2 = -3.85,
    a3 = 0.00024054,
    a4 = -116252,
    c1 = 284.6794,
    c2 = 471332,
    w_ref = 941530,
    z = -1,
    a_0 = 3.72,
    MM = 0.28704) "HUO3-";

  //752
  constant DataRecord UO2OH_aq(
    R = Modelica.Constants.R/UO2OH_aq.MM,
    G_ref = -1094534,
    H_ref = -1238046,
    S_ref = -58.158,
    a1 = 1.5644e-05,
    a2 = 563.46,
    a3 = 0.00021858,
    a4 = -118600,
    c1 = 85.4565,
    c2 = 84400,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.28603) "UO2OH,aq";

  //753
  constant DataRecord UO3_n1(
    R = Modelica.Constants.R/UO3_n1.MM,
    G_ref = -989934,
    H_ref = -1140558,
    S_ref = -81.588,
    a1 = 1.7317e-05,
    a2 = 973.7,
    a3 = 0.00020206,
    a4 = -120294,
    c1 = 116.3252,
    c2 = -69860,
    w_ref = 804290,
    z = -1,
    a_0 = 3.72,
    MM = 0.28603) "UO3-";

  //754
  constant DataRecord UO2OH_p1(
    R = Modelica.Constants.R/UO2OH_p1.MM,
    G_ref = -1160014,
    H_ref = -1261685,
    S_ref = 17.154,
    a1 = 1.9933e-05,
    a2 = 1612.05,
    a3 = 0.00017706,
    a4 = -122934,
    c1 = 55.3292,
    c2 = -90316,
    w_ref = 206060,
    z = 1,
    a_0 = 3.72,
    MM = 0.28704) "UO2OH+";

  //755
  constant DataRecord UO3_aq(
    R = Modelica.Constants.R/UO3_aq.MM,
    G_ref = -1130935,
    H_ref = -1253526,
    S_ref = -53.974,
    a1 = 1.5816e-05,
    a2 = 607.18,
    a3 = 0.00021646,
    a4 = -118780,
    c1 = 22.6827,
    c2 = -133783,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.23803) "UO3,aq";

  //756
  constant DataRecord HUO4_n1(
    R = Modelica.Constants.R/HUO4_n1.MM,
    G_ref = -1317123,
    H_ref = -1518374,
    S_ref = -84.517,
    a1 = 2.0654e-05,
    a2 = 1788.58,
    a3 = 0.00017003,
    a4 = -123662,
    c1 = 154.5365,
    c2 = 61388,
    w_ref = 809140,
    z = -1,
    a_0 = 3.72,
    MM = 0.30304) "HUO4-";

  //757
  constant DataRecord UO4_n2(
    R = Modelica.Constants.R/UO4_n2.MM,
    G_ref = -1238464,
    H_ref = -1448501,
    S_ref = -113.386,
    a1 = 2.1323e-05,
    a2 = 1950,
    a3 = 0.00016411,
    a4 = -124332,
    c1 = 56.5384,
    c2 = -505377,
    w_ref = 1515400,
    z = -2,
    a_0 = 3.72,
    MM = 0.30203) "UO4-2";

  //758
  constant DataRecord Fr_p1(
    R = Modelica.Constants.R/Fr_p1.MM,
    G_ref = -294972,
    H_ref = -261500,
    S_ref = 146.858,
    a1 = 1.9723e-05,
    a2 = 1559.96,
    a3 = 0.00017928,
    a4 = -122717,
    c1 = 10.6449,
    c2 = -182364,
    w_ref = 8490,
    z = 1,
    a_0 = 3.72,
    MM = 0.223) "Fr+";

  //759
  constant DataRecord OCN_n1(
    R = Modelica.Constants.R/OCN_n1.MM,
    G_ref = -97404,
    H_ref = -146022,
    S_ref = 106.692,
    a1 = 2.4096e-05,
    a2 = 2628.18,
    a3 = 0.00013719,
    a4 = -127135,
    c1 = 34.9402,
    c2 = -261626,
    w_ref = 519740,
    z = -1,
    a_0 = 3.72,
    MM = 0.042017) "OCN-";

  //760
  constant DataRecord SCN_n1(
    R = Modelica.Constants.R/SCN_n1.MM,
    G_ref = 92717,
    H_ref = 76442,
    S_ref = 144.348,
    a1 = 2.939e-05,
    a2 = 3919.86,
    a3 = 8.6642e-05,
    a4 = -132474,
    c1 = 44.942,
    c2 = -208782,
    w_ref = 463290,
    z = -1,
    a_0 = 3.72,
    MM = 0.058077) "SCN-";

  //761
  constant DataRecord SeCN_n1(
    R = Modelica.Constants.R/SeCN_n1.MM,
    G_ref = 81588,
    H_ref = 77530,
    S_ref = 194.974,
    a1 = 3.7134e-05,
    a2 = 5812.5,
    a3 = 1.1837e-05,
    a4 = -140298,
    c1 = 54.7589,
    c2 = -149975,
    w_ref = 386180,
    z = -1,
    a_0 = 3.72,
    MM = 0.10498) "SeCN-";

  //762
  constant DataRecord HNO2_aq(
    R = Modelica.Constants.R/HNO2_aq.MM,
    G_ref = -50626,
    H_ref = -119244,
    S_ref = 135.562,
    a1 = 2.4749e-05,
    a2 = 2786.13,
    a3 = 0.00013129,
    a4 = -127788,
    c1 = 36.4217,
    c2 = -69860,
    w_ref = -63050,
    z = 0,
    a_0 = 3.72,
    MM = 0.015015) "HNO2,aq";

  //763
  constant DataRecord H2N2O2_aq(
    R = Modelica.Constants.R/H2N2O2_aq.MM,
    G_ref = 35982,
    H_ref = -57321,
    S_ref = 217.568,
    a1 = 3.6123e-05,
    a2 = 5565.51,
    a3 = 2.1602e-05,
    a4 = -139277,
    c1 = 104.1833,
    c2 = 205422,
    w_ref = -187230,
    z = 0,
    a_0 = 3.72,
    MM = 0.030029) "H2N2O2,aq";

  //764
  constant DataRecord HN2O2_n1(
    R = Modelica.Constants.R/HN2O2_n1.MM,
    G_ref = 76149,
    H_ref = -39330,
    S_ref = 142.256,
    a1 = 2.7509e-05,
    a2 = 3462.39,
    a3 = 0.00010424,
    a4 = -130583,
    c1 = 91.5141,
    c2 = -47702,
    w_ref = 465760,
    z = -1,
    a_0 = 3.72,
    MM = 0.061021) "HN2O2-";

  //765
  constant DataRecord N2O2_n2(
    R = Modelica.Constants.R/N2O2_n2.MM,
    G_ref = 138909,
    H_ref = -10837,
    S_ref = 27.614,
    a1 = 1.385e-05,
    a2 = 124.93,
    a3 = 0.00023589,
    a4 = -116784,
    c1 = 0.4464,
    c2 = -632366,
    w_ref = 1303110,
    z = -2,
    a_0 = 3.72,
    MM = 0.060013) "N2O2-2";

  //766
  constant DataRecord N2H5_p1(
    R = Modelica.Constants.R/N2H5_p1.MM,
    G_ref = 82425,
    H_ref = -7531,
    S_ref = 151.042,
    a1 = 2.3309e-05,
    a2 = 2434.8,
    a3 = 0.00014507,
    a4 = -126336,
    c1 = 67.2134,
    c2 = 16217,
    w_ref = 2380,
    z = 1,
    a_0 = 3.72,
    MM = 0.033053) "N2H5+";

  //767
  constant DataRecord N2H6_p2(
    R = Modelica.Constants.R/N2H6_p2.MM,
    G_ref = 88282,
    H_ref = -16736,
    S_ref = 100.416,
    a1 = 1.8728e-05,
    a2 = 1316.62,
    a3 = 0.00018895,
    a4 = -121713,
    c1 = 91.033,
    c2 = 6841,
    w_ref = 290200,
    z = 2,
    a_0 = 3.72,
    MM = 0.034061) "N2H6+2";

  //768
  constant DataRecord H3PO2_aq(
    R = Modelica.Constants.R/H3PO2_aq.MM,
    G_ref = -523418,
    H_ref = -609190,
    S_ref = 154.39,
    a1 = 2.7343e-05,
    a2 = 3419.63,
    a3 = 0.0001064,
    a4 = -130407,
    c1 = 51.9402,
    c2 = -6791,
    w_ref = -91550,
    z = 0,
    a_0 = 3.72,
    MM = 0.033994) "H3PO2,aq";

  //769
  constant DataRecord H2PO2_n1(
    R = Modelica.Constants.R/H2PO2_n1.MM,
    G_ref = -512122,
    H_ref = -613793,
    S_ref = 100.834,
    a1 = 2.1768e-05,
    a2 = 2058.11,
    a3 = 0.00015992,
    a4 = -124779,
    c1 = 57.3472,
    c2 = -186623,
    w_ref = 528730,
    z = -1,
    a_0 = 3.72,
    MM = 0.064986) "H2PO2-";

  //770
  constant DataRecord H3PO3_aq(
    R = Modelica.Constants.R/H3PO3_aq.MM,
    G_ref = -856883,
    H_ref = -964830,
    S_ref = 182.422,
    a1 = 1.7867e-05,
    a2 = 1106.21,
    a3 = 0.00019726,
    a4 = -120842,
    c1 = 75.0028,
    c2 = 86956,
    w_ref = -134010,
    z = 0,
    a_0 = 3.72,
    MM = 0.033994) "H3PO3,aq";

  //771
  constant DataRecord H2PO3_n1(
    R = Modelica.Constants.R/H2PO3_n1.MM,
    G_ref = -846632,
    H_ref = -969433,
    S_ref = 132.633,
    a1 = 1.7426e-05,
    a2 = 999.93,
    a3 = 0.0002011,
    a4 = -120403,
    c1 = 83.5804,
    c2 = -80090,
    w_ref = 480780,
    z = -1,
    a_0 = 3.72,
    MM = 0.080986) "H2PO3-";

  //772
  constant DataRecord HPO3_n2(
    R = Modelica.Constants.R/HPO3_n2.MM,
    G_ref = -811696,
    H_ref = -969014,
    S_ref = 16.736,
    a1 = 1.4419e-05,
    a2 = 265.56,
    a3 = 0.00023001,
    a4 = -117365,
    c1 = 48.7532,
    c2 = -469579,
    w_ref = 1319090,
    z = -2,
    a_0 = 3.72,
    MM = 0.079978) "HPO3-2";

  //773
  constant DataRecord P2O7_n4(
    R = Modelica.Constants.R/P2O7_n4.MM,
    G_ref = -1919201,
    H_ref = -2271075,
    S_ref = -117.152,
    a1 = 2.9575e-05,
    a2 = 3965.3,
    a3 = 8.4822e-05,
    a4 = -132662,
    c1 = -49.0461,
    c2 = -1312487,
    w_ref = 2889850,
    z = -4,
    a_0 = 3.72,
    MM = 0.17394) "P2O7-4";

  //774
  constant DataRecord HP2O7_n3(
    R = Modelica.Constants.R/HP2O7_n3.MM,
    G_ref = -1972338,
    H_ref = -2274841,
    S_ref = 46.024,
    a1 = 3.4854e-05,
    a2 = 5253.35,
    a3 = 3.4342e-05,
    a4 = -137988,
    c1 = 134.552,
    c2 = -371568,
    w_ref = 1944300,
    z = -3,
    a_0 = 3.72,
    MM = 0.17495) "HP2O7-3";

  //775
  constant DataRecord H4P2O7_aq(
    R = Modelica.Constants.R/H4P2O7_aq.MM,
    G_ref = -2032169,
    H_ref = -2268565,
    S_ref = 267.776,
    a1 = 3.8901e-05,
    a2 = 6242.49,
    a3 = -4.728e-06,
    a4 = -142076,
    c1 = 145.7291,
    c2 = 374175,
    w_ref = -263260,
    z = 0,
    a_0 = 3.72,
    MM = 0.42682) "H4P2O7,aq";

  //776
  constant DataRecord AsO4_n3(
    R = Modelica.Constants.R/AsO4_n3.MM,
    G_ref = -646361,
    H_ref = -890209,
    S_ref = 176.314,
    a1 = 4.3129e-06,
    a2 = -2201.16,
    a3 = 0.00032673,
    a4 = -107169,
    c1 = -50.7737,
    c2 = -1116463,
    w_ref = 2258940,
    z = -3,
    a_0 = 3.72,
    MM = 0.13892) "AsO4-3";

  //777
  constant DataRecord H3AsO4_aq(
    R = Modelica.Constants.R/H3AsO4_aq.MM,
    G_ref = -766751,
    H_ref = -903451,
    S_ref = 183.05,
    a1 = 3.3203e-05,
    a2 = 4850.3,
    a3 = 5.0166e-05,
    a4 = -136319,
    c1 = 76.4856,
    c2 = 92922,
    w_ref = -136520,
    z = 0,
    a_0 = 3.72,
    MM = 0.077944) "H3AsO4,aq";

  //778
  constant DataRecord H2AsO4_1(
    R = Modelica.Constants.R/H2AsO4_1.MM,
    G_ref = -753651,
    H_ref = -911422,
    S_ref = 112.382,
    a1 = 0,
    a2 = 0,
    a3 = 0,
    a4 = 0,
    c1 = 0,
    c2 = 0,
    w_ref = 0,
    z = -1,
    a_0 = 3.72,
    MM = 0.14094) "H2AsO4-1";

  //779
  // constant DataRecord HAsO4_n2(
  //   R = Modelica.Constants.R/HAsO4_n2.MM,
  //   G_ref = -713732,
  //   H_ref = -908409,
  //   S_ref = 11.422,
  //   a1 = 0,
  //   a2 = 0,
  //   a3 = 0,
  //   a4 = 0,
  //   c1 = 0,
  //   c2 = 0,
  //   w_ref = 0,
  //   z = -2,
  //   a_0 = 3.72,
  //   MM = 0.13993) "HAsO4-2";

  //780
  constant DataRecord H3AsO3_aq(
    R = Modelica.Constants.R/H3AsO3_aq.MM,
    G_ref = -639805,
    H_ref = -742359,
    S_ref = 195.811,
    a1 = 0,
    a2 = 0,
    a3 = 0,
    a4 = 0,
    c1 = 0,
    c2 = 0,
    w_ref = 0,
    z = 0,
    a_0 = 3.72,
    MM = 0.077944) "H3AsO3,aq";

  //781
  constant DataRecord H2AsO3_1(
    R = Modelica.Constants.R/H2AsO3_1.MM,
    G_ref = -587660,
    H_ref = -714740,
    S_ref = 112.801,
    a1 = 0,
    a2 = 0,
    a3 = 0,
    a4 = 0,
    c1 = 0,
    c2 = 0,
    w_ref = 0,
    z = -1,
    a_0 = 3.72,
    MM = 0.12494) "H2AsO3-1";

  //782
  constant DataRecord HAsO3_n2(
    R = Modelica.Constants.R/HAsO3_n2.MM,
    G_ref = -507402,
    H_ref = 0,
    S_ref = 0,
    a1 = 0,
    a2 = 0,
    a3 = 0,
    a4 = 0,
    c1 = 0,
    c2 = 0,
    w_ref = 0,
    z = -2,
    a_0 = 3.72,
    MM = 0.12393) "HAsO3-2";

  //783
  constant DataRecord AsO3_n3(
    R = Modelica.Constants.R/AsO3_n3.MM,
    G_ref = -421802,
    H_ref = 0,
    S_ref = 0,
    a1 = 0,
    a2 = 0,
    a3 = 0,
    a4 = 0,
    c1 = 0,
    c2 = 0,
    w_ref = 0,
    z = -3,
    a_0 = 3.72,
    MM = 0.12292) "AsO3-3";

  //784
  constant DataRecord As3S4_HS2_n1(
    R = Modelica.Constants.R/As3S4_HS2_n1.MM,
    G_ref = -125599,
    H_ref = 0,
    S_ref = 0,
    a1 = 0,
    a2 = 0,
    a3 = 0,
    a4 = 0,
    c1 = 0,
    c2 = 0,
    w_ref = 0,
    z = -1,
    a_0 = 3.72,
    MM = 0.41914) "As3S4(HS)2-";

  //785
  constant DataRecord AsO2_n1(
    R = Modelica.Constants.R/AsO2_n1.MM,
    G_ref = -349992,
    H_ref = -429027,
    S_ref = 40.585,
    a1 = 1.3431e-05,
    a2 = 23.18,
    a3 = 0.00023976,
    a4 = -116365,
    c1 = 14.2691,
    c2 = -365602,
    w_ref = 620070,
    z = -1,
    a_0 = 3.72,
    MM = 0.10692) "AsO2-";

  //786
  constant DataRecord HAsO2_aq(
    R = Modelica.Constants.R/HAsO2_aq.MM,
    G_ref = -402668,
    H_ref = -456516,
    S_ref = 125.938,
    a1 = 2.3424e-05,
    a2 = 2464.92,
    a3 = 0.00014343,
    a4 = -126457,
    c1 = 28.447,
    c2 = -102249,
    w_ref = -48450,
    z = 0,
    a_0 = 3.72,
    MM = 0.16614) "HAsO2,aq";

  //787
  constant DataRecord HSbO2_aq(
    R = Modelica.Constants.R/HSbO2_aq.MM,
    G_ref = -407522,
    H_ref = -487854,
    S_ref = 46.442,
    a1 = 1.2437e-05,
    a2 = -217.69,
    a3 = 0.00024884,
    a4 = -115370,
    c1 = -37.4573,
    c2 = -369866,
    w_ref = 71920,
    z = 0,
    a_0 = 3.72,
    MM = 0.12175) "HSbO2,aq";

  //788
  constant DataRecord SbO2_n1(
    R = Modelica.Constants.R/SbO2_n1.MM,
    G_ref = -344762,
    H_ref = -460240,
    S_ref = -71.546,
    a1 = -2.1393e-06,
    a2 = -3776.85,
    a3 = 0.00038866,
    a4 = -100654,
    c1 = -9.5404,
    c2 = -502821,
    w_ref = 790150,
    z = -1,
    a_0 = 3.72,
    MM = 0.15375) "SbO2-";

  //789
  constant DataRecord Bi_p3(
    R = Modelica.Constants.R/Bi_p3.MM,
    G_ref = 95730,
    H_ref = 81002,
    S_ref = -188.28,
    a1 = -4.5844e-06,
    a2 = -4373.79,
    a3 = 0.00041215,
    a4 = -98186,
    c1 = 115.0713,
    c2 = -120144,
    w_ref = 947720,
    z = 3,
    a_0 = 3.72,
    MM = 0.20898) "Bi+3";

  //790
  constant DataRecord BiO_p1(
    R = Modelica.Constants.R/BiO_p1.MM,
    G_ref = -122591,
    H_ref = -126775,
    S_ref = 79.914,
    a1 = 1.2337e-05,
    a2 = -243.72,
    a3 = 0.00025025,
    a4 = -115261,
    c1 = -166.5897,
    c2 = -830947,
    w_ref = 110210,
    z = 1,
    a_0 = 3.72,
    MM = 0.22498) "BiO+";

  //791
  constant DataRecord BiOH_p2(
    R = Modelica.Constants.R/BiOH_p2.MM,
    G_ref = -135143,
    H_ref = -187443,
    S_ref = -81.588,
    a1 = 1.0332e-05,
    a2 = -733.83,
    a3 = 0.00026959,
    a4 = -113236,
    c1 = 41.2584,
    c2 = -254806,
    w_ref = 567020,
    z = 2,
    a_0 = 3.72,
    MM = 0.22599) "BiOH+2";

  //792
  constant DataRecord HBiO2_aq(
    R = Modelica.Constants.R/HBiO2_aq.MM,
    G_ref = -331791,
    H_ref = -361079,
    S_ref = 228.446,
    a1 = 2.1885e-05,
    a2 = 2087.48,
    a3 = 0.00015862,
    a4 = -124901,
    c1 = -339.4889,
    c2 = -1392603,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.24199) "HBiO2,aq";

  //793
  constant DataRecord BiO2_n1(
    R = Modelica.Constants.R/BiO2_n1.MM,
    G_ref = -258153,
    H_ref = -299574,
    S_ref = 188.698,
    a1 = 2.189e-05,
    a2 = 2088.65,
    a3 = 0.00015861,
    a4 = -124905,
    c1 = -250.1438,
    c2 = -1212770,
    w_ref = 395640,
    z = -1,
    a_0 = 3.72,
    MM = 0.24098) "BiO2-";

  //794
  constant DataRecord H2O2_aq(
    R = Modelica.Constants.R/H2O2_aq.MM,
    G_ref = -134014,
    H_ref = -191167,
    S_ref = 143.93,
    a1 = 2.5908e-05,
    a2 = 3069.8,
    a3 = 0.00012002,
    a4 = -128959,
    c1 = 43.3462,
    c2 = -41735,
    w_ref = -75690,
    z = 0,
    a_0 = 3.72,
    MM = 0.0020158) "H2O2,aq";

  //795
  constant DataRecord HClO2_aq(
    R = Modelica.Constants.R/HClO2_aq.MM,
    G_ref = 5858,
    H_ref = -51882,
    S_ref = 188.28,
    a1 = 3.2094e-05,
    a2 = 4579.6,
    a3 = 6.0781e-05,
    a4 = -135202,
    c1 = 80.0705,
    c2 = 107412,
    w_ref = -142880,
    z = 0,
    a_0 = 3.72,
    MM = 0.036458) "HClO2,aq";

  //796
  constant DataRecord HIO3_aq(
    R = Modelica.Constants.R/HIO3_aq.MM,
    G_ref = -132633,
    H_ref = -211292,
    S_ref = 166.942,
    a1 = 2.8596e-05,
    a2 = 3725.85,
    a3 = 9.4266e-05,
    a4 = -131670,
    c1 = 62.2052,
    c2 = 34966,
    w_ref = -110540,
    z = 0,
    a_0 = 3.72,
    MM = 0.12791) "HIO3,aq";

  //797
  constant DataRecord V_p2(
    R = Modelica.Constants.R/V_p2.MM,
    G_ref = -217568,
    H_ref = -225936,
    S_ref = -129.704,
    a1 = -6.7731e-06,
    a2 = -4907.37,
    a3 = 0.00043293,
    a4 = -95981,
    c1 = 57.9308,
    c2 = -219861,
    w_ref = 638850,
    z = 2,
    a_0 = 3.72,
    MM = 0.05094) "V+2";

  //798
  constant DataRecord V_p3(
    R = Modelica.Constants.R/V_p3.MM,
    G_ref = -242254,
    H_ref = -259408,
    S_ref = -230.12,
    a1 = -7.24e-06,
    a2 = -5020.01,
    a3 = 0.00043708,
    a4 = -95517,
    c1 = 98.1583,
    c2 = -198556,
    w_ref = 1008970,
    z = 3,
    a_0 = 3.72,
    MM = 0.05094) "V+3";

  //799
  constant DataRecord VO4_n3(
    R = Modelica.Constants.R/VO4_n3.MM,
    G_ref = -899142,
    H_ref = -1132190,
    S_ref = -147.277,
    a1 = 5.2631e-06,
    a2 = -1967.11,
    a3 = 0.00031714,
    a4 = -108135,
    c1 = 324.9106,
    c2 = 196899,
    w_ref = 2235260,
    z = -3,
    a_0 = 3.72,
    MM = 0.11494) "VO4-3";

  //800
  constant DataRecord Cr_p2(
    R = Modelica.Constants.R/Cr_p2.MM,
    G_ref = -164850,
    H_ref = -163176,
    S_ref = -101.253,
    a1 = -3.3623e-06,
    a2 = -4075.22,
    a3 = 0.00040036,
    a4 = -99420,
    c1 = 62.97,
    c2 = -189180,
    w_ref = 597770,
    z = 2,
    a_0 = 3.72,
    MM = 0.052) "Cr+2";

  //801
  constant DataRecord Cr_p3(
    R = Modelica.Constants.R/Cr_p3.MM,
    G_ref = -206271,
    H_ref = -251040,
    S_ref = -322.168,
    a1 = -1.1642e-05,
    a2 = -6096.46,
    a3 = 0.00047974,
    a4 = -91065,
    c1 = 75.276,
    c2 = -322135,
    w_ref = 1146540,
    z = 3,
    a_0 = 9,
    MM = 0.052) "Cr+3";

  //802
  constant DataRecord CrOH_p2(
    R = Modelica.Constants.R/CrOH_p2.MM,
    G_ref = -420492,
    H_ref = -496222,
    S_ref = -192.882,
    a1 = -1.108e-05,
    a2 = -5958.98,
    a3 = 0.00047424,
    a4 = -91634,
    c1 = 131.4307,
    c2 = 4284,
    w_ref = 736680,
    z = 2,
    a_0 = 3.72,
    MM = 0.069008) "CrOH+2";

  //803
  constant DataRecord CrO_p1(
    R = Modelica.Constants.R/CrO_p1.MM,
    G_ref = -388275,
    H_ref = -439320,
    S_ref = -110.039,
    a1 = -4.907e-06,
    a2 = -4451.15,
    a3 = 0.00041489,
    a4 = -97868,
    c1 = -13.3474,
    c2 = -389468,
    w_ref = 394840,
    z = 1,
    a_0 = 3.72,
    MM = 0.068) "CrO+";

  //804
  constant DataRecord HCrO2_aq(
    R = Modelica.Constants.R/HCrO2_aq.MM,
    G_ref = -577810,
    H_ref = -658143,
    S_ref = 25.104,
    a1 = 6.4262e-06,
    a2 = -1684.44,
    a3 = 0.00030632,
    a4 = -109307,
    c1 = -67.7988,
    c2 = -448274,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.085008) "HCrO2,aq";

  //805
  constant DataRecord CrO2_n1(
    R = Modelica.Constants.R/CrO2_n1.MM,
    G_ref = -524255,
    H_ref = -620487,
    S_ref = -28.033,
    a1 = 5.0262e-06,
    a2 = -2025.18,
    a3 = 0.00031945,
    a4 = -107897,
    c1 = 41.8488,
    c2 = -303386,
    w_ref = 725130,
    z = -1,
    a_0 = 3.72,
    MM = 0.084) "CrO2-";

  //806
  constant DataRecord Zr_p4(
    R = Modelica.Constants.R/Zr_p4.MM,
    G_ref = -557602,
    H_ref = -628981,
    S_ref = -461.495,
    a1 = -1.8271e-05,
    a2 = -7714.33,
    a3 = 0.00054313,
    a4 = -84379,
    c1 = 166.2818,
    c2 = -153385,
    w_ref = 1607370,
    z = 4,
    a_0 = 3.72,
    MM = 0.09122) "Zr+4";

  //807
  constant DataRecord ZrOH_p3(
    R = Modelica.Constants.R/ZrOH_p3.MM,
    G_ref = -796634,
    H_ref = -889100,
    S_ref = -299.574,
    a1 = 3.8823e-06,
    a2 = -2305.09,
    a3 = 0.0003306,
    a4 = -106738,
    c1 = 237.4152,
    c2 = 251446,
    w_ref = 1115200,
    z = 3,
    a_0 = 3.72,
    MM = 0.10823) "ZrOH+3";

  //808
  constant DataRecord ZrO_p2(
    R = Modelica.Constants.R/ZrO_p2.MM,
    G_ref = -784918,
    H_ref = -854791,
    S_ref = -223.426,
    a1 = 7.9555e-06,
    a2 = -1312.27,
    a3 = 0.00029191,
    a4 = -110843,
    c1 = 97.7851,
    c2 = -126110,
    w_ref = 778680,
    z = 2,
    a_0 = 3.72,
    MM = 0.10722) "ZrO+2";

  //809
  constant DataRecord HZrO2_p1(
    R = Modelica.Constants.R/HZrO2_p1.MM,
    G_ref = -1002905,
    H_ref = -1110015,
    S_ref = -207.945,
    a1 = 8.8308e-06,
    a2 = -1100.35,
    a3 = 0.000284,
    a4 = -111721,
    c1 = 295.1155,
    c2 = 634119,
    w_ref = 546430,
    z = 1,
    a_0 = 3.72,
    MM = 0.12423) "HZrO2+";

  //810
  constant DataRecord ZrO2_aq(
    R = Modelica.Constants.R/ZrO2_aq.MM,
    G_ref = -976546,
    H_ref = -1103739,
    S_ref = -182.841,
    a1 = 7.8002e-06,
    a2 = -1349.93,
    a3 = 0.00029338,
    a4 = -110688,
    c1 = 160.7351,
    c2 = 346050,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.12322) "ZrO2,aq";

  //811
  constant DataRecord HZrO3_n1(
    R = Modelica.Constants.R/HZrO3_n1.MM,
    G_ref = -1177796,
    H_ref = -1394527,
    S_ref = -136.817,
    a1 = 1.4905e-05,
    a2 = 383.59,
    a3 = 0.00022552,
    a4 = -117855,
    c1 = 231.5915,
    c2 = 304286,
    w_ref = 886970,
    z = -1,
    a_0 = 3.72,
    MM = 0.14023) "HZrO3-";

  //812
  constant DataRecord HNbO3_aq(
    R = Modelica.Constants.R/HNbO3_aq.MM,
    G_ref = -990771,
    H_ref = -1079263,
    S_ref = 130.122,
    a1 = 2.3895e-05,
    a2 = 2579.1,
    a3 = 0.00013913,
    a4 = -126930,
    c1 = 84.6419,
    c2 = 140649,
    w_ref = -197070,
    z = 0,
    a_0 = 3.72,
    MM = 0.14192) "HNbO3,aq";

  //813
  constant DataRecord NbO3_n1(
    R = Modelica.Constants.R/NbO3_n1.MM,
    G_ref = -950186,
    H_ref = -1068175,
    S_ref = 13.389,
    a1 = 2.1303e-05,
    a2 = 1944.93,
    a3 = 0.00016432,
    a4 = -124311,
    c1 = -13.9645,
    c2 = -477252,
    w_ref = 662290,
    z = -1,
    a_0 = 3.72,
    MM = 0.14091) "NbO3-";

  //814
  constant DataRecord TcO4_n1(
    R = Modelica.Constants.R/TcO4_n1.MM,
    G_ref = -632202,
    H_ref = -724669,
    S_ref = 198.74,
    a1 = 3.6039e-05,
    a2 = 7834.79,
    a3 = -0.00055874,
    a4 = -148658,
    c1 = 79.2868,
    c2 = -248948,
    w_ref = 380790,
    z = -1,
    a_0 = 3.72,
    MM = 0.16291) "TcO4-";

  //815
  constant DataRecord Pm_p2(
    R = Modelica.Constants.R/Pm_p2.MM,
    G_ref = -387857,
    H_ref = -369866,
    S_ref = -7.95,
    a1 = 2.176e-07,
    a2 = -3201.22,
    a3 = 0.00036608,
    a4 = -103035,
    c1 = 39.4564,
    c2 = -224978,
    w_ref = 454300,
    z = 2,
    a_0 = 3.72,
    MM = 0.145) "Pm+2";

  //816
  constant DataRecord Pm_p3(
    R = Modelica.Constants.R/Pm_p3.MM,
    G_ref = -661072,
    H_ref = -686176,
    S_ref = -209.2,
    a1 = -1.3586e-05,
    a2 = -6571.35,
    a3 = 0.00049842,
    a4 = -89102,
    c1 = 6.5174,
    c2 = -507080,
    w_ref = 977760,
    z = 3,
    a_0 = 3.72,
    MM = 0.145) "Pm+3";

  //817
  constant DataRecord Pm_p4(
    R = Modelica.Constants.R/Pm_p4.MM,
    G_ref = -143093,
    H_ref = -210874,
    S_ref = -417.145,
    a1 = -1.7923e-05,
    a2 = -7628.23,
    a3 = 0.00053955,
    a4 = -84734,
    c1 = 168.2776,
    c2 = -125261,
    w_ref = 1541180,
    z = 4,
    a_0 = 3.72,
    MM = 0.145) "Pm+4";

  //818
  constant DataRecord Hf_p4(
    R = Modelica.Constants.R/Hf_p4.MM,
    G_ref = -554798,
    H_ref = -628646,
    S_ref = -465.261,
    a1 = -1.831e-05,
    a2 = -7725.25,
    a3 = 0.0005439,
    a4 = -84333,
    c1 = 166.0688,
    c2 = -155942,
    w_ref = 1613020,
    z = 4,
    a_0 = 3.72,
    MM = 0.17849) "Hf+4";

  //819
  constant DataRecord HfOH_p3(
    R = Modelica.Constants.R/HfOH_p3.MM,
    G_ref = -790776,
    H_ref = -886171,
    S_ref = -304.595,
    a1 = 3.7853e-06,
    a2 = -2329.23,
    a3 = 0.00033165,
    a4 = -106642,
    c1 = 241.3222,
    c2 = 263379,
    w_ref = 1120350,
    z = 3,
    a_0 = 3.72,
    MM = 0.1955) "HfOH+3";

  //820
  constant DataRecord HfO_p2(
    R = Modelica.Constants.R/HfO_p2.MM,
    G_ref = -778224,
    H_ref = -840566,
    S_ref = -193.301,
    a1 = 8.4429e-06,
    a2 = -1194.57,
    a3 = 0.00028762,
    a4 = -111332,
    c1 = 73.8062,
    c2 = -195999,
    w_ref = 736680,
    z = 2,
    a_0 = 3.72,
    MM = 0.10722) "HfO+2";

  //821
  constant DataRecord HHfO2_p1(
    R = Modelica.Constants.R/HHfO2_p1.MM,
    G_ref = -994955,
    H_ref = -1105413,
    S_ref = -213.384,
    a1 = 8.6312e-06,
    a2 = -1147.29,
    a3 = 0.00028544,
    a4 = -111525,
    c1 = 303.2705,
    c2 = 659687,
    w_ref = 555130,
    z = 1,
    a_0 = 3.72,
    MM = 0.2115) "HHfO2+";

  //822
  constant DataRecord HfO2_aq(
    R = Modelica.Constants.R/HfO2_aq.MM,
    G_ref = -968178,
    H_ref = -1098718,
    S_ref = -188.698,
    a1 = 7.5714e-06,
    a2 = -1407.46,
    a3 = 0.00029597,
    a4 = -110449,
    c1 = 167.8466,
    c2 = 370765,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.21049) "HfO2,aq";

  //823
  constant DataRecord HHfO3_n1(
    R = Modelica.Constants.R/HHfO3_n1.MM,
    G_ref = -1168173,
    H_ref = -1388251,
    S_ref = -143.511,
    a1 = 1.4605e-05,
    a2 = 310.75,
    a3 = 0.00022828,
    a4 = -117554,
    c1 = 241.603,
    c2 = 334967,
    w_ref = 899810,
    z = -1,
    a_0 = 3.72,
    MM = 0.2275) "HHfO3-";

  //824
  constant DataRecord La_p2(
    R = Modelica.Constants.R/La_p2.MM,
    G_ref = -325097,
    H_ref = -302085,
    S_ref = 2.929,
    a1 = 4.452e-07,
    a2 = -3143.36,
    a3 = 0.00036332,
    a4 = -103274,
    c1 = 36.6313,
    c2 = -229237,
    w_ref = 436980,
    z = 2,
    a_0 = 3.72,
    MM = 0.13891) "La+2";

  //825
  constant DataRecord LaCO3_p1(
    R = Modelica.Constants.R/LaCO3_p1.MM,
    G_ref = -1254782,
    H_ref = -1400845,
    S_ref = -184.933,
    a1 = -3.4464e-06,
    a2 = -4094.96,
    a3 = 0.00040099,
    a4 = -99341,
    c1 = -84.0599,
    c2 = -673277,
    w_ref = 513590,
    z = 1,
    a_0 = 3.72,
    MM = 0.19892) "LaCO3+";

  //826
  constant DataRecord LaHCO3_p2(
    R = Modelica.Constants.R/LaHCO3_p2.MM,
    G_ref = -1284488,
    H_ref = -1392854,
    S_ref = -57.321,
    a1 = 3.7606e-06,
    a2 = -2335.84,
    a3 = 0.00033203,
    a4 = -106613,
    c1 = 159.0744,
    c2 = 166201,
    w_ref = 531080,
    z = 2,
    a_0 = 3.72,
    MM = 0.19993) "LaHCO3+2";

  //827
  constant DataRecord LaOH_p2(
    R = Modelica.Constants.R/LaOH_p2.MM,
    G_ref = -874038,
    H_ref = -910355,
    S_ref = -27.614,
    a1 = 1.1199e-05,
    a2 = -522.41,
    a3 = 0.00026133,
    a4 = -114110,
    c1 = -2.1192,
    c2 = -379238,
    w_ref = 484760,
    z = 2,
    a_0 = 3.72,
    MM = 0.15592) "LaOH+2";

  //828
  constant DataRecord LaO_p1(
    R = Modelica.Constants.R/LaO_p1.MM,
    G_ref = -819646,
    H_ref = -836256,
    S_ref = 38.493,
    a1 = 1.1629e-05,
    a2 = -415.81,
    a3 = 0.00025679,
    a4 = -114550,
    c1 = -132.9826,
    c2 = -733786,
    w_ref = 171540,
    z = 1,
    a_0 = 3.72,
    MM = 0.15491) "LaO+";

  //829
  constant DataRecord LaO2H_aq(
    R = Modelica.Constants.R/LaO2H_aq.MM,
    G_ref = -1001231,
    H_ref = -1043908,
    S_ref = 184.096,
    a1 = 2.0339e-05,
    a2 = 1709.29,
    a3 = 0.00017365,
    a4 = -123336,
    c1 = -279.9042,
    c2 = -1185495,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17091) "LaO2H,aq";

  //830
  constant DataRecord LaO2_n1(
    R = Modelica.Constants.R/LaO2_n1.MM,
    G_ref = -927593,
    H_ref = -983031,
    S_ref = 141.419,
    a1 = 2.0473e-05,
    a2 = 1742.09,
    a3 = 0.00017229,
    a4 = -123470,
    c1 = -186.152,
    c2 = -1013336,
    w_ref = 467440,
    z = -1,
    a_0 = 3.72,
    MM = 0.17091) "LaO2-";

  //831
  constant DataRecord LaCl_p2(
    R = Modelica.Constants.R/LaCl_p2.MM,
    G_ref = -819227,
    H_ref = -862322,
    S_ref = -107.529,
    a1 = 1.142e-07,
    a2 = -3226.37,
    a3 = 0.00036704,
    a4 = -102931,
    c1 = 70.4707,
    c2 = -165661,
    w_ref = 605720,
    z = 2,
    a_0 = 3.72,
    MM = 0.17436) "LaCl+2";

  //832
  constant DataRecord LaCl2_p1(
    R = Modelica.Constants.R/LaCl2_p1.MM,
    G_ref = -948513,
    H_ref = -1024662,
    S_ref = -40.166,
    a1 = 1.3299e-05,
    a2 = -7.41,
    a3 = 0.00024064,
    a4 = -116240,
    c1 = 28.5533,
    c2 = -211221,
    w_ref = 293010,
    z = 1,
    a_0 = 3.72,
    MM = 0.20981) "LaCl2+";

  //833
  constant DataRecord LaCl3_aq(
    R = Modelica.Constants.R/LaCl3_aq.MM,
    G_ref = -1077798,
    H_ref = -1198632,
    S_ref = -13.807,
    a1 = 2.8134e-05,
    a2 = 3614.39,
    a3 = 9.8395e-05,
    a4 = -131210,
    c1 = -105.8954,
    c2 = -580685,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13891) "LaCl3,aq";

  //834
  constant DataRecord LaCl4_n1(
    R = Modelica.Constants.R/LaCl4_n1.MM,
    G_ref = -1206666,
    H_ref = -1385741,
    S_ref = -32.635,
    a1 = 4.8332e-05,
    a2 = 8546.66,
    a3 = -9.5558e-05,
    a4 = -151603,
    c1 = -236.7495,
    c2 = -1274061,
    w_ref = 732450,
    z = -1,
    a_0 = 3.72,
    MM = 0.28071) "LaCl4-";

  //835
  constant DataRecord LaNO3_p2(
    R = Modelica.Constants.R/LaNO3_p2.MM,
    G_ref = -800399,
    H_ref = -945584,
    S_ref = -158.155,
    a1 = 6.7944e-06,
    a2 = -1595.53,
    a3 = 0.00030305,
    a4 = -109675,
    c1 = 136.9687,
    c2 = 40606,
    w_ref = 683330,
    z = 2,
    a_0 = 3.72,
    MM = 0.20092) "LaNO3+2";

  //836
  constant DataRecord LaF_p2(
    R = Modelica.Constants.R/LaF_p2.MM,
    G_ref = -989934,
    H_ref = -1018386,
    S_ref = -68.199,
    a1 = -1.1031e-05,
    a2 = -5945.59,
    a3 = 0.00047339,
    a4 = -91688,
    c1 = 74.9467,
    c2 = -130704,
    w_ref = 545130,
    z = 2,
    a_0 = 3.72,
    MM = 0.15791) "LaF+2";

  //837
  constant DataRecord LaF2_p1(
    R = Modelica.Constants.R/LaF2_p1.MM,
    G_ref = -1287835,
    H_ref = -1360637,
    S_ref = -50.208,
    a1 = -9.7922e-06,
    a2 = -5643.67,
    a3 = 0.00046167,
    a4 = -92939,
    c1 = 53.6644,
    c2 = -127742,
    w_ref = 304890,
    z = 1,
    a_0 = 3.72,
    MM = 0.17691) "LaF2+";

  //838
  constant DataRecord LaF3_aq(
    R = Modelica.Constants.R/LaF3_aq.MM,
    G_ref = -1581134,
    H_ref = -1716277,
    S_ref = -93.303,
    a1 = -8.5797e-06,
    a2 = -5347.82,
    a3 = 0.0004501,
    a4 = -94161,
    c1 = -64.016,
    c2 = -435128,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13891) "LaF3,aq";

  //839
  constant DataRecord LaF4_n1(
    R = Modelica.Constants.R/LaF4_n1.MM,
    G_ref = -1872340,
    H_ref = -2092418,
    S_ref = -210.874,
    a1 = -2.6204e-06,
    a2 = -3893.63,
    a3 = 0.00039312,
    a4 = -100173,
    c1 = -148.7127,
    c2 = -1052858,
    w_ref = 997210,
    z = -1,
    a_0 = 3.72,
    MM = 0.21491) "LaF4-";

  //840
  constant DataRecord LaH2PO4_p2(
    R = Modelica.Constants.R/LaH2PO4_p2.MM,
    G_ref = -1830918,
    H_ref = -2020035,
    S_ref = -126.357,
    a1 = 7.9467e-06,
    a2 = -1315.2,
    a3 = 0.00029223,
    a4 = -110834,
    c1 = 175.259,
    c2 = 189297,
    w_ref = 634630,
    z = 2,
    a_0 = 3.72,
    MM = 0.2359) "LaH2PO4+2";

  //841
  constant DataRecord LaSO4_p1(
    R = Modelica.Constants.R/LaSO4_p1.MM,
    G_ref = -1451430,
    H_ref = -1600798,
    S_ref = -66.944,
    a1 = 6.7266e-06,
    a2 = -1612.1,
    a3 = 0.00030369,
    a4 = -109608,
    c1 = -113.5429,
    c2 = 629813,
    w_ref = 330330,
    z = 1,
    a_0 = 3.72,
    MM = 0.23497) "LaSO4+";

  //842
  constant DataRecord CeCO3_p1(
    R = Modelica.Constants.R/CeCO3_p1.MM,
    G_ref = -1246359,
    H_ref = -1391682,
    S_ref = -169.87,
    a1 = -3.8325e-06,
    a2 = -4189.61,
    a3 = 0.00040474,
    a4 = -98947,
    c1 = -89.5983,
    c2 = -685209,
    w_ref = 490740,
    z = 1,
    a_0 = 3.72,
    MM = 0.20013) "CeCO3+";

  //843
  constant DataRecord CeHCO3_p2(
    R = Modelica.Constants.R/CeHCO3_p2.MM,
    G_ref = -1274028,
    H_ref = -1381557,
    S_ref = -39.748,
    a1 = 2.9041e-06,
    a2 = -2544.83,
    a3 = 0.0003402,
    a4 = -105751,
    c1 = 148.8634,
    c2 = 139356,
    w_ref = 504090,
    z = 2,
    a_0 = 3.72,
    MM = 0.20114) "CeHCO3+2";

  //844
  constant DataRecord CeOH_p2(
    R = Modelica.Constants.R/CeOH_p2.MM,
    G_ref = -865251,
    H_ref = -900439,
    S_ref = -11.297,
    a1 = 1.1517e-05,
    a2 = -443.76,
    a3 = 0.00025804,
    a4 = -114437,
    c1 = -15.4139,
    c2 = -417592,
    w_ref = 460240,
    z = 2,
    a_0 = 3.72,
    MM = 0.15713) "CeOH+2";

  //845
  constant DataRecord CeO_p1(
    R = Modelica.Constants.R/CeO_p1.MM,
    G_ref = -819646,
    H_ref = -834750,
    S_ref = 56.066,
    a1 = 1.1886e-05,
    a2 = -352.71,
    a3 = 0.00025427,
    a4 = -114813,
    c1 = -147.3446,
    c2 = -775550,
    w_ref = 146060,
    z = 1,
    a_0 = 3.72,
    MM = 0.15612) "CeO+";

  //846
  constant DataRecord CeO2H_aq(
    R = Modelica.Constants.R/CeO2H_aq.MM,
    G_ref = -1001231,
    H_ref = -1043908,
    S_ref = 202.924,
    a1 = 2.0969e-05,
    a2 = 1865.31,
    a3 = 0.00016701,
    a4 = -123980,
    c1 = -305.4056,
    c2 = -1274133,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17212) "CeO2H,aq";

  //847
  constant DataRecord CeO2_n1(
    R = Modelica.Constants.R/CeO2_n1.MM,
    G_ref = -929266,
    H_ref = -982487,
    S_ref = 161.502,
    a1 = 2.1115e-05,
    a2 = 1900.29,
    a3 = 0.0001658,
    a4 = -124127,
    c1 = -213.7053,
    c2 = -1099417,
    w_ref = 437190,
    z = -1,
    a_0 = 3.72,
    MM = 0.17212) "CeO2-";

  //848
  constant DataRecord CeCl_p2(
    R = Modelica.Constants.R/CeCl_p2.MM,
    G_ref = -809186,
    H_ref = -852699,
    S_ref = -92.466,
    a1 = -7.305e-07,
    a2 = -3430.59,
    a3 = 0.00037457,
    a4 = -102085,
    c1 = 60.576,
    c2 = -192510,
    w_ref = 582160,
    z = 2,
    a_0 = 3.72,
    MM = 0.17557) "CeCl+2";

  //849
  constant DataRecord CeCl2_p1(
    R = Modelica.Constants.R/CeCl2_p1.MM,
    G_ref = -938890,
    H_ref = -1013783,
    S_ref = -21.338,
    a1 = 1.2347e-05,
    a2 = -241.12,
    a3 = 0.00025011,
    a4 = -115273,
    c1 = 10.7918,
    c2 = -263600,
    w_ref = 263800,
    z = 1,
    a_0 = 3.72,
    MM = 0.21102) "CeCl2+";

  //850
  constant DataRecord CeCl3_aq(
    R = Modelica.Constants.R/CeCl3_aq.MM,
    G_ref = -1067757,
    H_ref = -1186164,
    S_ref = 10.878,
    a1 = 2.7184e-05,
    a2 = 3382.97,
    a3 = 0.00010736,
    a4 = -130256,
    c1 = -131.3671,
    c2 = -669222,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14012) "CeCl3,aq";

  //851
  constant DataRecord CeCl4_n1(
    R = Modelica.Constants.R/CeCl4_n1.MM,
    G_ref = -1197042,
    H_ref = -1370678,
    S_ref = 1.255,
    a1 = 4.7094e-05,
    a2 = 8244.4,
    a3 = -8.3693e-05,
    a4 = -150352,
    c1 = -280.564,
    c2 = -1409368,
    w_ref = 679400,
    z = -1,
    a_0 = 3.72,
    MM = 0.28192) "CeCl4-";

  //852
  constant DataRecord CeH2PO4_p2(
    R = Modelica.Constants.R/CeH2PO4_p2.MM,
    G_ref = -1820458,
    H_ref = -2008738,
    S_ref = -108.366,
    a1 = 7.0835e-06,
    a2 = -1525.07,
    a3 = 0.00030029,
    a4 = -109964,
    c1 = 164.8739,
    c2 = 162448,
    w_ref = 605720,
    z = 2,
    a_0 = 3.72,
    MM = 0.23711) "CeH2PO4+2";

  //853
  constant DataRecord CeNO3_p2(
    R = Modelica.Constants.R/CeNO3_p2.MM,
    G_ref = -794542,
    H_ref = -933869,
    S_ref = -134.725,
    a1 = 5.9078e-06,
    a2 = -1809.71,
    a3 = 0.00031095,
    a4 = -108788,
    c1 = 125.9409,
    c2 = 13757,
    w_ref = 647470,
    z = 2,
    a_0 = 3.72,
    MM = 0.20213) "CeNO3+2";

  //854
  constant DataRecord CeF_p2(
    R = Modelica.Constants.R/CeF_p2.MM,
    G_ref = -981985,
    H_ref = -1012528,
    S_ref = -59.831,
    a1 = -1.1833e-05,
    a2 = -6142.95,
    a3 = 0.00048151,
    a4 = -90872,
    c1 = 66.2469,
    c2 = -157553,
    w_ref = 534550,
    z = 2,
    a_0 = 3.72,
    MM = 0.15912) "CeF+2";

  //855
  constant DataRecord CeF2_p1(
    R = Modelica.Constants.R/CeF2_p1.MM,
    G_ref = -1281141,
    H_ref = -1356034,
    S_ref = -42.258,
    a1 = -1.0685e-05,
    a2 = -5861.45,
    a3 = 0.00047017,
    a4 = -92040,
    c1 = 37.5004,
    c2 = -180125,
    w_ref = 293010,
    z = 1,
    a_0 = 3.72,
    MM = 0.17812) "CeF2+";

  //856
  constant DataRecord CeF3_aq(
    R = Modelica.Constants.R/CeF3_aq.MM,
    G_ref = -1575694,
    H_ref = -1712511,
    S_ref = -82.425,
    a1 = -9.5299e-06,
    a2 = -5579.28,
    a3 = 0.00045906,
    a4 = -93203,
    c1 = -89.4886,
    c2 = -523661,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14012) "CeF3,aq";

  //857
  constant DataRecord CeF4_n1(
    R = Modelica.Constants.R/CeF4_n1.MM,
    G_ref = -1868156,
    H_ref = -2087398,
    S_ref = -192.882,
    a1 = -3.7635e-06,
    a2 = -4171.16,
    a3 = 0.00040371,
    a4 = -99027,
    c1 = -189.9364,
    c2 = -1188164,
    w_ref = 972320,
    z = -1,
    a_0 = 3.72,
    MM = 0.21612) "CeF4-";

  //858
  constant DataRecord CeBr_p2(
    R = Modelica.Constants.R/CeBr_p2.MM,
    G_ref = -782358,
    H_ref = -818846,
    S_ref = -104.642,
    a1 = 3.3773e-06,
    a2 = -2428.6,
    a3 = 0.00033547,
    a4 = -106228,
    c1 = 57.9183,
    c2 = -208012,
    w_ref = 601740,
    z = 2,
    a_0 = 3.72,
    MM = 0.22003) "CeBr+2";

  //859
  constant DataRecord CeIO3_p2(
    R = Modelica.Constants.R/CeIO3_p2.MM,
    G_ref = -815010,
    H_ref = -942898,
    S_ref = -121.21,
    a1 = 4.061e-06,
    a2 = -2260.49,
    a3 = 0.00032865,
    a4 = -106926,
    c1 = 124.8903,
    c2 = 16920,
    w_ref = 626220,
    z = 2,
    a_0 = 3.72,
    MM = 0.31502) "CeIO3+2";

  //860
  constant DataRecord CeClO4_p2(
    R = Modelica.Constants.R/CeClO4_p2.MM,
    G_ref = -695573,
    H_ref = -878749,
    S_ref = -150.875,
    a1 = 1.4685e-05,
    a2 = 328.99,
    a3 = 0.00022784,
    a4 = -117629,
    c1 = 177.5878,
    c2 = 186175,
    w_ref = 669650,
    z = 2,
    a_0 = 3.72,
    MM = 0.23957) "CeClO4+2";

  //861
  constant DataRecord CeSO4_p1(
    R = Modelica.Constants.R/CeSO4_p1.MM,
    G_ref = -1399548,
    H_ref = -1590757,
    S_ref = -52.718,
    a1 = 6.368e-06,
    a2 = -1699.62,
    a3 = 0.00030709,
    a4 = -109244,
    c1 = -93.8408,
    c2 = -641742,
    w_ref = 308950,
    z = 1,
    a_0 = 3.72,
    MM = 0.23618) "CeSO4+";

  //862
  constant DataRecord PrCO3_p1(
    R = Modelica.Constants.R/PrCO3_p1.MM,
    G_ref = -1251434,
    H_ref = -1399590,
    S_ref = -174.891,
    a1 = -4.4028e-06,
    a2 = -4327.39,
    a3 = 0.00040985,
    a4 = -98378,
    c1 = -93.5701,
    c2 = -701402,
    w_ref = 498190,
    z = 1,
    a_0 = 3.72,
    MM = 0.20092) "PrCO3+";

  //863
  constant DataRecord PrHCO3_p2(
    R = Modelica.Constants.R/PrHCO3_p2.MM,
    G_ref = -1278212,
    H_ref = -1409171,
    S_ref = -118.407,
    a1 = 1.8355e-06,
    a2 = -2805.33,
    a3 = 0.00035035,
    a4 = -104671,
    c1 = 149.25,
    c2 = 102918,
    w_ref = 622040,
    z = 2,
    a_0 = 3.72,
    MM = 0.20193) "PrHCO3+2";

  //864
  constant DataRecord PrCl_p2(
    R = Modelica.Constants.R/PrCl_p2.MM,
    G_ref = -813370,
    H_ref = -858975,
    S_ref = -97.487,
    a1 = -2.1715e-06,
    a2 = -3784.09,
    a3 = 0.00038882,
    a4 = -100625,
    c1 = 50.8059,
    c2 = -228944,
    w_ref = 589900,
    z = 2,
    a_0 = 3.72,
    MM = 0.17636) "PrCl+2";

  //865
  constant DataRecord PrCl2_p1(
    R = Modelica.Constants.R/PrCl2_p1.MM,
    G_ref = -943074,
    H_ref = -1020059,
    S_ref = -27.614,
    a1 = 1.0748e-05,
    a2 = -631.53,
    a3 = 0.00026544,
    a4 = -113658,
    c1 = -8.6877,
    c2 = -334691,
    w_ref = 274340,
    z = 1,
    a_0 = 3.72,
    MM = 0.21181) "PrCl2+";

  //866
  constant DataRecord PrCl3_aq(
    R = Modelica.Constants.R/PrCl3_aq.MM,
    G_ref = -1072359,
    H_ref = -1193277,
    S_ref = 2.51,
    a1 = 2.5363e-05,
    a2 = 2936.04,
    a3 = 0.00012541,
    a4 = -128407,
    c1 = -165.937,
    c2 = -789374,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14091) "PrCl3,aq";

  //867
  constant DataRecord PrCl4_n1(
    R = Modelica.Constants.R/PrCl4_n1.MM,
    G_ref = -1201226,
    H_ref = -1378628,
    S_ref = -10.042,
    a1 = 4.5127e-05,
    a2 = 7762.41,
    a3 = -6.4358e-05,
    a4 = -148360,
    c1 = -331.6895,
    c2 = -1592999,
    w_ref = 697930,
    z = -1,
    a_0 = 3.72,
    MM = 0.28271) "PrCl4-";

  //868
  constant DataRecord PrH2PO4_p2(
    R = Modelica.Constants.R/PrH2PO4_p2.MM,
    G_ref = -1824224,
    H_ref = -2014596,
    S_ref = -114.223,
    a1 = 5.6438e-06,
    a2 = -1876.27,
    a3 = 0.000314,
    a4 = -108512,
    c1 = 155.1369,
    c2 = 126014,
    w_ref = 613830,
    z = 2,
    a_0 = 3.72,
    MM = 0.2379) "PrH2PO4+2";

  //869
  constant DataRecord PrNO3_p2(
    R = Modelica.Constants.R/PrNO3_p2.MM,
    G_ref = -794960,
    H_ref = -940982,
    S_ref = -142.256,
    a1 = 4.4702e-06,
    a2 = -2161.45,
    a3 = 0.00032492,
    a4 = -107332,
    c1 = 116.2637,
    c2 = -22673,
    w_ref = 656220,
    z = 2,
    a_0 = 3.72,
    MM = 0.20292) "PrNO3+2";

  //870
  constant DataRecord PrF_p2(
    R = Modelica.Constants.R/PrF_p2.MM,
    G_ref = -986169,
    H_ref = -1018386,
    S_ref = -62.76,
    a1 = -1.3288e-05,
    a2 = -6497.96,
    a3 = 0.00049539,
    a4 = -89408,
    c1 = 56.0874,
    c2 = -193987,
    w_ref = 538060,
    z = 2,
    a_0 = 3.72,
    MM = 0.15991) "PrF+2";

  //871
  constant DataRecord PrF2_p1(
    R = Modelica.Constants.R/PrF2_p1.MM,
    G_ref = -1285743,
    H_ref = -1362310,
    S_ref = -44.769,
    a1 = -1.2307e-05,
    a2 = -6258.93,
    a3 = 0.00048613,
    a4 = -90395,
    c1 = 17.4063,
    c2 = -251212,
    w_ref = 296900,
    z = 1,
    a_0 = 3.72,
    MM = 0.17891) "PrF2+";

  //872
  constant DataRecord PrF3_aq(
    R = Modelica.Constants.R/PrF3_aq.MM,
    G_ref = -1580715,
    H_ref = -1718787,
    S_ref = -85.772,
    a1 = -1.1351e-05,
    a2 = -6026.3,
    a3 = 0.00047712,
    a4 = -91358,
    c1 = -124.0577,
    c2 = -643817,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14091) "PrF3,aq";

  //873
  constant DataRecord PrF4_n1(
    R = Modelica.Constants.R/PrF4_n1.MM,
    G_ref = -1873177,
    H_ref = -2094929,
    S_ref = -199.158,
    a1 = -5.7656e-06,
    a2 = -4662.23,
    a3 = 0.00042348,
    a4 = -96993,
    c1 = -242.021,
    c2 = -1371796,
    w_ref = 980440,
    z = -1,
    a_0 = 3.72,
    MM = 0.21691) "PrF4-";

  //874
  constant DataRecord PrOH_p2(
    R = Modelica.Constants.R/PrOH_p2.MM,
    G_ref = -870272,
    H_ref = -908388,
    S_ref = -16.736,
    a1 = 1.1433e-05,
    a2 = -463.63,
    a3 = 0.00025866,
    a4 = -114353,
    c1 = -10.9027,
    c2 = -404806,
    w_ref = 469280,
    z = 2,
    a_0 = 3.72,
    MM = 0.15792) "PrOH+2";

  //875
  constant DataRecord PrO_p1(
    R = Modelica.Constants.R/PrO_p1.MM,
    G_ref = -818809,
    H_ref = -836967,
    S_ref = 50.208,
    a1 = 1.1799e-05,
    a2 = -374.26,
    a3 = 0.00025517,
    a4 = -114721,
    c1 = -142.6736,
    c2 = -761911,
    w_ref = 154180,
    z = 1,
    a_0 = 3.72,
    MM = 0.15691) "PrO+";

  //876
  constant DataRecord PrO2H_aq(
    R = Modelica.Constants.R/PrO2H_aq.MM,
    G_ref = -1002905,
    H_ref = -1046418,
    S_ref = 196.648,
    a1 = 2.0797e-05,
    a2 = 1821.59,
    a3 = 0.00016913,
    a4 = -123800,
    c1 = -309.621,
    c2 = -1245154,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17291) "PrO2H,aq";

  //877
  constant DataRecord PrO2_n1(
    R = Modelica.Constants.R/PrO2_n1.MM,
    G_ref = -940145,
    H_ref = -996671,
    S_ref = 154.808,
    a1 = 2.092e-05,
    a2 = 1852.05,
    a3 = 0.00016782,
    a4 = -123926,
    c1 = -204.4332,
    c2 = -1070439,
    w_ref = 447310,
    z = -1,
    a_0 = 3.72,
    MM = 0.17291) "PrO2-";

  //878
  constant DataRecord PrSO4_p1(
    R = Modelica.Constants.R/PrSO4_p1.MM,
    G_ref = -1403732,
    H_ref = -1596196,
    S_ref = -54.81,
    a1 = -5.8095e-06,
    a2 = -1836.02,
    a3 = 0.00031245,
    a4 = -108679,
    c1 = -98.1181,
    c2 = -657938,
    w_ref = 313090,
    z = 1,
    a_0 = 3.72,
    MM = 0.23697) "PrSO4+";

  //879
  constant DataRecord NdCO3_p1(
    R = Modelica.Constants.R/NdCO3_p1.MM,
    G_ref = -1243903,
    H_ref = -1390469,
    S_ref = -172.381,
    a1 = -4.6857e-06,
    a2 = -4396.88,
    a3 = 0.00041269,
    a4 = -98094,
    c1 = -100.8783,
    c2 = -724414,
    w_ref = 490740,
    z = 1,
    a_0 = 3.72,
    MM = 0.20425) "NdCO3+";

  //880
  constant DataRecord NdHCO3_p2(
    R = Modelica.Constants.R/NdHCO3_p2.MM,
    G_ref = -1269426,
    H_ref = -1377373,
    S_ref = -42.677,
    a1 = 8.1e-07,
    a2 = -3056.2,
    a3 = 0.00036031,
    a4 = -103633,
    c1 = 123.7874,
    c2 = 51145,
    w_ref = 507350,
    z = 2,
    a_0 = 3.72,
    MM = 0.20526) "NdHCO3+2";

  //881
  constant DataRecord NdOH_p2(
    R = Modelica.Constants.R/NdOH_p2.MM,
    G_ref = -862741,
    H_ref = -899142,
    S_ref = -13.807,
    a1 = 1.147e-05,
    a2 = -454.97,
    a3 = 0.00025838,
    a4 = -114386,
    c1 = -13.421,
    c2 = -411626,
    w_ref = 463250,
    z = 2,
    a_0 = 3.72,
    MM = 0.16125) "NdOH+2";

  //882
  constant DataRecord NdO_p1(
    R = Modelica.Constants.R/NdO_p1.MM,
    G_ref = -811696,
    H_ref = -828139,
    S_ref = 53.137,
    a1 = 1.1842e-05,
    a2 = -364.68,
    a3 = 0.00025504,
    a4 = -114763,
    c1 = -145.0128,
    c2 = -768731,
    w_ref = 150080,
    z = 1,
    a_0 = 3.72,
    MM = 0.16024) "NdO+";

  //883
  constant DataRecord NdO2H_aq(
    R = Modelica.Constants.R/NdO2H_aq.MM,
    G_ref = -995792,
    H_ref = -1037632,
    S_ref = 199.995,
    a1 = 2.0854e-05,
    a2 = 1835.31,
    a3 = 0.00016863,
    a4 = -123855,
    c1 = -301.2371,
    c2 = -1259643,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17624) "NdO2H,aq";

  //884
  constant DataRecord NdO2_n1(
    R = Modelica.Constants.R/NdO2_n1.MM,
    G_ref = -934706,
    H_ref = -989390,
    S_ref = 158.155,
    a1 = 2.1016e-05,
    a2 = 1876.52,
    a3 = 0.00016663,
    a4 = -124026,
    c1 = -209.1113,
    c2 = -1084928,
    w_ref = 441790,
    z = -1,
    a_0 = 3.72,
    MM = 0.17624) "NdO2-";

  //885
  constant DataRecord NdCl_p2(
    R = Modelica.Constants.R/NdCl_p2.MM,
    G_ref = -805002,
    H_ref = -849352,
    S_ref = -94.977,
    a1 = -2.8225e-06,
    a2 = -3942.5,
    a3 = 0.00039493,
    a4 = -99972,
    c1 = 35.5523,
    c2 = -280721,
    w_ref = 586010,
    z = 2,
    a_0 = 3.72,
    MM = 0.17969) "NdCl+2";

  //886
  constant DataRecord NdCl2_p1(
    R = Modelica.Constants.R/NdCl2_p1.MM,
    G_ref = -934706,
    H_ref = -1010436,
    S_ref = -24.686,
    a1 = 1.0014e-05,
    a2 = -809.77,
    a3 = 0.0002722,
    a4 = -112922,
    c1 = -38.4049,
    c2 = -435713,
    w_ref = 267270,
    z = 1,
    a_0 = 3.72,
    MM = 0.21514) "NdCl2+";

  //887
  constant DataRecord NdCl3_aq(
    R = Modelica.Constants.R/NdCl3_aq.MM,
    G_ref = -1063991,
    H_ref = -1182942,
    S_ref = 6.694,
    a1 = 2.4571e-05,
    a2 = 2743.95,
    a3 = 0.00013266,
    a4 = -127612,
    c1 = -215.0614,
    c2 = -960123,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14424) "NdCl3,aq";

  //888
  constant DataRecord NdCl4_n1(
    R = Modelica.Constants.R/NdCl4_n1.MM,
    G_ref = -1192858,
    H_ref = -1368168,
    S_ref = -4.602,
    a1 = 4.4213e-05,
    a2 = 7539.78,
    a3 = -5.5748e-05,
    a4 = -147440,
    c1 = -407.6379,
    c2 = -1853951,
    w_ref = 688480,
    z = -1,
    a_0 = 3.72,
    MM = 0.28604) "NdCl4-";

  //889
  constant DataRecord NdH2PO4_p2(
    R = Modelica.Constants.R/NdH2PO4_p2.MM,
    G_ref = -1815400,
    H_ref = -2004454,
    S_ref = -111.169,
    a1 = 4.9919e-06,
    a2 = -2034.47,
    a3 = 0.00032004,
    a4 = -107859,
    c1 = 139.8665,
    c2 = 74237,
    w_ref = 609780,
    z = 2,
    a_0 = 3.72,
    MM = 0.24123) "NdH2PO4+2";

  //890
  constant DataRecord NdNO3_p2(
    R = Modelica.Constants.R/NdNO3_p2.MM,
    G_ref = -787366,
    H_ref = -931300,
    S_ref = -138.449,
    a1 = 3.8175e-06,
    a2 = -2322.04,
    a3 = 0.0003315,
    a4 = -106671,
    c1 = 100.9624,
    c2 = -74450,
    w_ref = 651830,
    z = 2,
    a_0 = 3.72,
    MM = 0.20625) "NdNO3+2";

  //891
  constant DataRecord NdF_p2(
    R = Modelica.Constants.R/NdF_p2.MM,
    G_ref = -978638,
    H_ref = -1009181,
    S_ref = -61.086,
    a1 = -1.3938e-05,
    a2 = -6656.66,
    a3 = 0.00050164,
    a4 = -88751,
    c1 = 40.8672,
    c2 = -245760,
    w_ref = 534550,
    z = 2,
    a_0 = 3.72,
    MM = 0.16324) "NdF+2";

  //892
  constant DataRecord NdF2_p1(
    R = Modelica.Constants.R/NdF2_p1.MM,
    G_ref = -1278630,
    H_ref = -1353524,
    S_ref = -43.514,
    a1 = -1.3018e-05,
    a2 = -6432.98,
    a3 = 0.00049308,
    a4 = -89676,
    c1 = -11.6575,
    c2 = -352234,
    w_ref = 296900,
    z = 1,
    a_0 = 3.72,
    MM = 0.18224) "NdF2+";

  //893
  constant DataRecord NdF3_aq(
    R = Modelica.Constants.R/NdF3_aq.MM,
    G_ref = -1573602,
    H_ref = -1710838,
    S_ref = -84.098,
    a1 = -1.2143e-05,
    a2 = -6218.26,
    a3 = 0.00048437,
    a4 = -90563,
    c1 = -173.1829,
    c2 = -814562,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14424) "NdF3,aq";

  //894
  constant DataRecord NdF4_n1(
    R = Modelica.Constants.R/NdF4_n1.MM,
    G_ref = -1866482,
    H_ref = -2086561,
    S_ref = -195.811,
    a1 = -6.648e-06,
    a2 = -4877.46,
    a3 = 0.00043191,
    a4 = -96106,
    c1 = -317.0987,
    c2 = -1632747,
    w_ref = 980440,
    z = -1,
    a_0 = 3.72,
    MM = 0.22024) "NdF4-";

  //895
  constant DataRecord NdSO4_p1(
    R = Modelica.Constants.R/NdSO4_p1.MM,
    G_ref = -1437204,
    H_ref = -1586154,
    S_ref = -51.463,
    a1 = 5.5091e-06,
    a2 = -1909.33,
    a3 = 0.00031534,
    a4 = -108378,
    c1 = -105.1205,
    c2 = -680950,
    w_ref = 313130,
    z = 1,
    a_0 = 3.72,
    MM = 0.2403) "NdSO4+";

  //896
  constant DataRecord SmCO3_p1(
    R = Modelica.Constants.R/SmCO3_p1.MM,
    G_ref = -1238464,
    H_ref = -1386327,
    S_ref = -178.238,
    a1 = -4.3744e-06,
    a2 = -4321.78,
    a3 = 0.00040995,
    a4 = -98403,
    c1 = -100.4365,
    c2 = -725263,
    w_ref = 498190,
    z = 1,
    a_0 = 3.72,
    MM = 0.21036) "SmCO3+";

  //897
  constant DataRecord SmOH_p2(
    R = Modelica.Constants.R/SmOH_p2.MM,
    G_ref = -857302,
    H_ref = -895250,
    S_ref = -20.502,
    a1 = 1.1329e-05,
    a2 = -488.52,
    a3 = 0.00025952,
    a4 = -114248,
    c1 = -8.1684,
    c2 = -396283,
    w_ref = 472330,
    z = 2,
    a_0 = 3.72,
    MM = 0.16736) "SmOH+2";

  //898
  constant DataRecord SmO_p1(
    R = Modelica.Constants.R/SmO_p1.MM,
    G_ref = -808767,
    H_ref = -826884,
    S_ref = 46.024,
    a1 = 1.1763e-05,
    a2 = -383.13,
    a3 = 0.00025554,
    a4 = -114683,
    c1 = -139.1477,
    c2 = -751685,
    w_ref = 160540,
    z = 1,
    a_0 = 3.72,
    MM = 0.16635) "SmO+";

  //899
  constant DataRecord SmO2H_aq(
    R = Modelica.Constants.R/SmO2H_aq.MM,
    G_ref = -992026,
    H_ref = -1036377,
    S_ref = 192.046,
    a1 = 2.0625e-05,
    a2 = 1780.38,
    a3 = 0.00017057,
    a4 = -123629,
    c1 = -290.9382,
    c2 = -1223849,
    w_ref = -125520,
    z = 0,
    a_0 = 3.72,
    MM = 0.18235) "SmO2H,aq";

  //900
  constant DataRecord SmO2_n1(
    R = Modelica.Constants.R/SmO2_n1.MM,
    G_ref = -940145,
    H_ref = -996712,
    S_ref = 150.206,
    a1 = 2.077e-05,
    a2 = 1815.56,
    a3 = 0.00016927,
    a4 = -123775,
    c1 = -198.1881,
    c2 = -1050837,
    w_ref = 453880,
    z = -1,
    a_0 = 3.72,
    MM = 0.18235) "SmO2-";

  //901
  constant DataRecord SmHCO3_p2(
    R = Modelica.Constants.R/SmHCO3_p2.MM,
    G_ref = -1262731,
    H_ref = -1371934,
    S_ref = -49.79,
    a1 = 1.5456e-06,
    a2 = -2875.54,
    a3 = 0.00035298,
    a4 = -104382,
    c1 = 124.1589,
    c2 = 49225,
    w_ref = 517390,
    z = 2,
    a_0 = 3.72,
    MM = 0.21137) "SmHCO3+2";

  //902
  constant DataRecord SmCl_p2(
    R = Modelica.Constants.R/SmCl_p2.MM,
    G_ref = -798726,
    H_ref = -843913,
    S_ref = -100.834,
    a1 = -2.0945e-06,
    a2 = -3765.1,
    a3 = 0.00038804,
    a4 = -100705,
    c1 = 35.7184,
    c2 = -282638,
    w_ref = 593790,
    z = 2,
    a_0 = 3.72,
    MM = 0.1858) "SmCl+2";

  //903
  constant DataRecord SmCl2_p1(
    R = Modelica.Constants.R/SmCl2_p1.MM,
    G_ref = -928011,
    H_ref = -1005415,
    S_ref = -32.217,
    a1 = 1.0831e-05,
    a2 = -611.58,
    a3 = 0.00026475,
    a4 = -113742,
    c1 = -38.4949,
    c2 = -439454,
    w_ref = 277980,
    z = 1,
    a_0 = 3.72,
    MM = 0.22125) "SmCl2+";

  //904
  constant DataRecord SmCl3_aq(
    R = Modelica.Constants.R/SmCl3_aq.MM,
    G_ref = -1057297,
    H_ref = -1178633,
    S_ref = -2.929,
    a1 = 2.5442e-05,
    a2 = 2956.96,
    a3 = 0.00012423,
    a4 = -128495,
    c1 = -216.8814,
    c2 = -966445,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15035) "SmCl3,aq";

  //905
  constant DataRecord SmCl4_n1(
    R = Modelica.Constants.R/SmCl4_n1.MM,
    G_ref = -1186164,
    H_ref = -1364821,
    S_ref = -17.991,
    a1 = 4.5249e-05,
    a2 = 7793.16,
    a3 = -6.5823e-05,
    a4 = -148486,
    c1 = -408.6408,
    c2 = -1863616,
    w_ref = 707810,
    z = -1,
    a_0 = 3.72,
    MM = 0.29215) "SmCl4-";

  //906
  constant DataRecord SmH2PO4_p2(
    R = Modelica.Constants.R/SmH2PO4_p2.MM,
    G_ref = -1808743,
    H_ref = -1999115,
    S_ref = -118.407,
    a1 = 5.7354e-06,
    a2 = -1853.3,
    a3 = 0.00031297,
    a4 = -108608,
    c1 = 140.4464,
    c2 = 72320,
    w_ref = 622040,
    z = 2,
    a_0 = 3.72,
    MM = 0.24734) "SmH2PO4+2";

  //907
  constant DataRecord SmNO3_p2(
    R = Modelica.Constants.R/SmNO3_p2.MM,
    G_ref = -781153,
    H_ref = -927174,
    S_ref = -148.114,
    a1 = 4.5639e-06,
    a2 = -2139.03,
    a3 = 0.00032417,
    a4 = -107428,
    c1 = 101.6344,
    c2 = -76366,
    w_ref = 665130,
    z = 2,
    a_0 = 3.72,
    MM = 0.21236) "SmNO3+2";

  //908
  constant DataRecord SmF_p2(
    R = Modelica.Constants.R/SmF_p2.MM,
    G_ref = -972362,
    H_ref = -1003742,
    S_ref = -64.434,
    a1 = -1.3212e-05,
    a2 = -6478.63,
    a3 = 0.00049445,
    a4 = -89487,
    c1 = 40.9647,
    c2 = -247680,
    w_ref = 541580,
    z = 2,
    a_0 = 3.72,
    MM = 0.16935) "SmF+2";

  //909
  constant DataRecord SmF2_p1(
    R = Modelica.Constants.R/SmF2_p1.MM,
    G_ref = -1272773,
    H_ref = -1348085,
    S_ref = -46.861,
    a1 = -1.2222e-05,
    a2 = -6236.71,
    a3 = 0.0004849,
    a4 = -90487,
    c1 = -12.3696,
    c2 = -355975,
    w_ref = 300870,
    z = 1,
    a_0 = 3.72,
    MM = 0.18835) "SmF2+";

  //910
  constant DataRecord SmF3_aq(
    R = Modelica.Constants.R/SmF3_aq.MM,
    G_ref = -1568163,
    H_ref = -1705817,
    S_ref = -88.701,
    a1 = -1.1272e-05,
    a2 = -6005.25,
    a3 = 0.00047594,
    a4 = -91446,
    c1 = -175.0029,
    c2 = -820884,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15035) "SmF3,aq";

  //911
  constant DataRecord SmF4_n1(
    R = Modelica.Constants.R/SmF4_n1.MM,
    G_ref = -1861043,
    H_ref = -2082377,
    S_ref = -202.924,
    a1 = -5.6492e-06,
    a2 = -4632.65,
    a3 = 0.00042207,
    a4 = -97119,
    c1 = -319.1149,
    c2 = -1642412,
    w_ref = 988760,
    z = -1,
    a_0 = 3.72,
    MM = 0.22635) "SmF4-";

  //912
  constant DataRecord SmSO4_p1(
    R = Modelica.Constants.R/SmSO4_p1.MM,
    G_ref = -1430928,
    H_ref = -1580715,
    S_ref = -55.647,
    a1 = -5.8095e-06,
    a2 = -1836.02,
    a3 = 0.00031245,
    a4 = -108679,
    c1 = -104.9841,
    c2 = -681800,
    w_ref = 313090,
    z = 1,
    a_0 = 3.72,
    MM = 0.24641) "SmSO4+";

  //913
  constant DataRecord EuCO3_p1(
    R = Modelica.Constants.R/EuCO3_p1.MM,
    G_ref = -1147671,
    H_ref = -1302354,
    S_ref = -189.954,
    a1 = -4.1179e-06,
    a2 = -4258.43,
    a3 = 0.00040728,
    a4 = -98663,
    c1 = -81.8566,
    c2 = -668164,
    w_ref = 521540,
    z = 1,
    a_0 = 3.72,
    MM = 0.21197) "EuCO3+";

  //914
  constant DataRecord EuOH_p2(
    R = Modelica.Constants.R/EuOH_p2.MM,
    G_ref = -766509,
    H_ref = -811696,
    S_ref = -33.472,
    a1 = 1.1116e-05,
    a2 = -542.62,
    a3 = 0.00026217,
    a4 = -114027,
    c1 = 2.4384,
    c2 = -366456,
    w_ref = 494340,
    z = 2,
    a_0 = 3.72,
    MM = 0.16897) "EuOH+2";

  //915
  constant DataRecord EuO_p1(
    R = Modelica.Constants.R/EuO_p1.MM,
    G_ref = -718393,
    H_ref = -743957,
    S_ref = 32.217,
    a1 = 1.1488e-05,
    a2 = -449.49,
    a3 = 0.000258,
    a4 = -114411,
    c1 = -128.204,
    c2 = -720150,
    w_ref = 180830,
    z = 1,
    a_0 = 3.72,
    MM = 0.16796) "EuO+";

  //916
  constant DataRecord EuO2H_aq(
    R = Modelica.Constants.R/EuO2H_aq.MM,
    G_ref = -903744,
    H_ref = -954789,
    S_ref = 177.402,
    a1 = 2.011e-05,
    a2 = 1654.35,
    a3 = 0.00017559,
    a4 = -123110,
    c1 = -271.5667,
    c2 = -1156520,
    w_ref = -125520,
    z = 0,
    a_0 = 3.72,
    MM = 0.18396) "EuO2H,aq";

  //917
  constant DataRecord EuO2_n1(
    R = Modelica.Constants.R/EuO2_n1.MM,
    G_ref = -851862,
    H_ref = -916547,
    S_ref = 134.306,
    a1 = 2.0279e-05,
    a2 = 1696.24,
    a3 = 0.00017384,
    a4 = -123282,
    c1 = -177.0874,
    c2 = -985211,
    w_ref = 477980,
    z = -1,
    a_0 = 3.72,
    MM = 0.18396) "EuO2-";

  //918
  constant DataRecord EuHCO3_p2(
    R = Modelica.Constants.R/EuHCO3_p2.MM,
    G_ref = -1170683,
    H_ref = -1286580,
    S_ref = -63.597,
    a1 = 2.0619e-06,
    a2 = -2749.72,
    a3 = 0.0003481,
    a4 = -104901,
    c1 = 163.0287,
    c2 = 177707,
    w_ref = 538060,
    z = 2,
    a_0 = 3.72,
    MM = 0.21298) "EuHCO3+2";

  //919
  constant DataRecord EuCl_p2(
    R = Modelica.Constants.R/EuCl_p2.MM,
    G_ref = -707514,
    H_ref = -758559,
    S_ref = -112.55,
    a1 = -1.5803e-06,
    a2 = -3638.74,
    a3 = 0.0003829,
    a4 = -101228,
    c1 = 74.5267,
    c2 = -154155,
    w_ref = 613830,
    z = 2,
    a_0 = 3.72,
    MM = 0.18741) "EuCl+2";

  //920
  constant DataRecord EuCl2_p1(
    R = Modelica.Constants.R/EuCl2_p1.MM,
    G_ref = -836800,
    H_ref = -920898,
    S_ref = -46.442,
    a1 = 1.1406e-05,
    a2 = -469.44,
    a3 = 0.00025886,
    a4 = -114328,
    c1 = 35.7351,
    c2 = -188770,
    w_ref = 300870,
    z = 1,
    a_0 = 3.72,
    MM = 0.22286) "EuCl2+";

  //921
  constant DataRecord EuCl3_aq(
    R = Modelica.Constants.R/EuCl3_aq.MM,
    G_ref = -965667,
    H_ref = -1095371,
    S_ref = -22.175,
    a1 = 2.5996e-05,
    a2 = 3091.18,
    a3 = 0.00011922,
    a4 = -129047,
    c1 = -94.9781,
    c2 = -542744,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15196) "EuCl3,aq";

  //922
  constant DataRecord EuCl4_n1(
    R = Modelica.Constants.R/EuCl4_n1.MM,
    G_ref = -1094534,
    H_ref = -1283651,
    S_ref = -43.932,
    a1 = 4.6001e-05,
    a2 = 7977.21,
    a3 = -7.3107e-05,
    a4 = -149247,
    c1 = -218.6609,
    c2 = -1216071,
    w_ref = 747680,
    z = -1,
    a_0 = 3.72,
    MM = 0.29376) "EuCl4-";

  //923
  constant DataRecord EuF_p1(
    R = Modelica.Constants.R/EuF_p1.MM,
    G_ref = -814206,
    H_ref = -846005,
    S_ref = 7.113,
    a1 = 1.1089e-05,
    a2 = -548.52,
    a3 = 0.00026221,
    a4 = -114001,
    c1 = 172.2473,
    c2 = 311624,
    w_ref = 219910,
    z = 1,
    a_0 = 3.72,
    MM = 0.17096) "EuF+";

  //924
  constant DataRecord EuF2_aq(
    R = Modelica.Constants.R/EuF2_aq.MM,
    G_ref = -1092024,
    H_ref = -1180725,
    S_ref = -17.154,
    a1 = 1.5006e-05,
    a2 = 409.36,
    a3 = 0.00022425,
    a4 = -117960,
    c1 = 267.2944,
    c2 = 716422,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15196) "EuF2,aq";

  //925
  constant DataRecord EuF3_n1(
    R = Modelica.Constants.R/EuF3_n1.MM,
    G_ref = -1371097,
    H_ref = -1530089,
    S_ref = -85.772,
    a1 = 2.3023e-05,
    a2 = 2366.3,
    a3 = 0.00014746,
    a4 = -126051,
    c1 = 464.447,
    c2 = 1138554,
    w_ref = 809140,
    z = -1,
    a_0 = 3.72,
    MM = 0.20896) "EuF3-";

  //926
  constant DataRecord EuF4_n2(
    R = Modelica.Constants.R/EuF4_n2.MM,
    G_ref = -1651006,
    H_ref = -1906649,
    S_ref = -243.09,
    a1 = 3.1919e-05,
    a2 = 4536.75,
    a3 = 6.2492e-05,
    a4 = -135026,
    c1 = 674.2311,
    c2 = 1578029,
    w_ref = 1713810,
    z = -2,
    a_0 = 3.72,
    MM = 0.22796) "EuF4-2";

  //927
  constant DataRecord EuCl_p1(
    R = Modelica.Constants.R/EuCl_p1.MM,
    G_ref = -673624,
    H_ref = -686176,
    S_ref = 81.588,
    a1 = 2.1649e-05,
    a2 = 2029.2,
    a3 = 0.00016103,
    a4 = -124658,
    c1 = 151.7859,
    c2 = 276667,
    w_ref = 106980,
    z = 1,
    a_0 = 3.72,
    MM = 0.18741) "EuCl+";

  //928
  constant DataRecord EuCl2_aq(
    R = Modelica.Constants.R/EuCl2_aq.MM,
    G_ref = -810022,
    H_ref = -856046,
    S_ref = 146.44,
    a1 = 3.8138e-05,
    a2 = 6055.92,
    a3 = 2.682e-06,
    a4 = -141306,
    c1 = 243.2774,
    c2 = 632943,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15196) "EuCl2,aq";

  //929
  constant DataRecord EuCl3_n1(
    R = Modelica.Constants.R/EuCl3_n1.MM,
    G_ref = -945584,
    H_ref = -1032611,
    S_ref = 187.025,
    a1 = 5.835e-05,
    a2 = 10992.25,
    a3 = -0.00019159,
    a4 = -161712,
    c1 = 384.7485,
    c2 = 992993,
    w_ref = 398610,
    z = -1,
    a_0 = 3.72,
    MM = 0.25831) "EuCl3-";

  //930
  constant DataRecord EuCl4_n2(
    R = Modelica.Constants.R/EuCl4_n2.MM,
    G_ref = -1081564,
    H_ref = -1215870,
    S_ref = 203.342,
    a1 = 8.1475e-05,
    a2 = 16637.93,
    a3 = -0.00041331,
    a4 = -185050,
    c1 = 548.1228,
    c2 = 1356825,
    w_ref = 1035750,
    z = -2,
    a_0 = 3.72,
    MM = 0.29376) "EuCl4-2";

  //931
  constant DataRecord EuH2PO4_p2(
    R = Modelica.Constants.R/EuH2PO4_p2.MM,
    G_ref = -1717532,
    H_ref = -1914598,
    S_ref = -132.633,
    a1 = 6.2534e-06,
    a2 = -1725.31,
    a3 = 0.0003076,
    a4 = -109135,
    c1 = 179.3568,
    c2 = 200803,
    w_ref = 643160,
    z = 2,
    a_0 = 3.72,
    MM = 0.24895) "EuH2PO4+2";

  //932
  constant DataRecord EuNO3_p2(
    R = Modelica.Constants.R/EuNO3_p2.MM,
    G_ref = -690360,
    H_ref = -844331,
    S_ref = -166.105,
    a1 = 5.1036e-06,
    a2 = -2006.27,
    a3 = 0.00031873,
    a4 = -107976,
    c1 = 141.1389,
    c2 = 52112,
    w_ref = 692700,
    z = 2,
    a_0 = 3.72,
    MM = 0.21397) "EuNO3+2";

  //933
  constant DataRecord EuF_p2(
    R = Modelica.Constants.R/EuF_p2.MM,
    G_ref = -881569,
    H_ref = -917133,
    S_ref = -71.128,
    a1 = -1.2741e-05,
    a2 = -6364.37,
    a3 = 0.0004901,
    a4 = -89960,
    c1 = 78.5872,
    c2 = -119198,
    w_ref = 548730,
    z = 2,
    a_0 = 3.72,
    MM = 0.17096) "EuF+2";

  //934
  constant DataRecord EuF2_p1(
    R = Modelica.Constants.R/EuF2_p1.MM,
    G_ref = -1181980,
    H_ref = -1262313,
    S_ref = -53.137,
    a1 = -1.1698e-05,
    a2 = -6110.48,
    a3 = 0.00048038,
    a4 = -91010,
    c1 = 60.4994,
    c2 = -105295,
    w_ref = 308950,
    z = 1,
    a_0 = 3.72,
    MM = 0.18996) "EuF2+";

  //935
  constant DataRecord EuF3_aq(
    R = Modelica.Constants.R/EuF3_aq.MM,
    G_ref = -1477789,
    H_ref = -1620463,
    S_ref = -96.65,
    a1 = -1.0718e-05,
    a2 = -5871.16,
    a3 = 0.00047093,
    a4 = -91998,
    c1 = -53.0996,
    c2 = -397183,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15196) "EuF3,aq";

  //936
  constant DataRecord EuF4_n1(
    R = Modelica.Constants.R/EuF4_n1.MM,
    G_ref = -1770669,
    H_ref = -1999115,
    S_ref = -216.731,
    a1 = -4.9735e-06,
    a2 = -4467.68,
    a3 = 0.00041558,
    a4 = -97801,
    c1 = -131.2303,
    c2 = -994867,
    w_ref = 1005920,
    z = -1,
    a_0 = 3.72,
    MM = 0.22796) "EuF4-";

  //937
  constant DataRecord EuSO4_p1(
    R = Modelica.Constants.R/EuSO4_p1.MM,
    G_ref = -1339717,
    H_ref = -1452685,
    S_ref = -64.015,
    a1 = -6.0245e-06,
    a2 = -1783.51,
    a3 = 0.00031039,
    a4 = -108897,
    c1 = -87.3741,
    c2 = -624696,
    w_ref = 325930,
    z = 1,
    a_0 = 3.72,
    MM = 0.24802) "EuSO4+";

  //938
  constant DataRecord GdCO3_p1(
    R = Modelica.Constants.R/GdCO3_p1.MM,
    G_ref = -1236372,
    H_ref = -1381640,
    S_ref = -170.707,
    a1 = -3.9874e-06,
    a2 = -4227.35,
    a3 = 0.00040625,
    a4 = -98793,
    c1 = -82.9779,
    c2 = -662197,
    w_ref = 490740,
    z = 1,
    a_0 = 3.72,
    MM = 0.21726) "GdCO3+";

  //939
  constant DataRecord GdOH_p2(
    R = Modelica.Constants.R/GdOH_p2.MM,
    G_ref = -855628,
    H_ref = -890774,
    S_ref = -12.134,
    a1 = 1.146e-05,
    a2 = -457.56,
    a3 = 0.00025851,
    a4 = -114378,
    c1 = -14.6783,
    c2 = -415036,
    w_ref = 460240,
    z = 2,
    a_0 = 3.72,
    MM = 0.17426) "GdOH+2";

  //940
  constant DataRecord GdO_p1(
    R = Modelica.Constants.R/GdO_p1.MM,
    G_ref = -807512,
    H_ref = -822700,
    S_ref = 54.81,
    a1 = 1.1893e-05,
    a2 = -351.83,
    a3 = 0.00025439,
    a4 = -114813,
    c1 = -146.4245,
    c2 = -772994,
    w_ref = 148070,
    z = 1,
    a_0 = 3.72,
    MM = 0.17325) "GdO+";

  //941
  constant DataRecord GdO2H_aq(
    R = Modelica.Constants.R/GdO2H_aq.MM,
    G_ref = -993700,
    H_ref = -1034285,
    S_ref = 201.669,
    a1 = 2.0969e-05,
    a2 = 1865.31,
    a3 = 0.00016701,
    a4 = -123980,
    c1 = -303.6894,
    c2 = -1268166,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.18925) "GdO2H,aq";

  //942
  constant DataRecord GdO2_n1(
    R = Modelica.Constants.R/GdO2_n1.MM,
    G_ref = -941400,
    H_ref = -994662,
    S_ref = 160.247,
    a1 = 2.1064e-05,
    a2 = 1887.44,
    a3 = 0.00016639,
    a4 = -124072,
    c1 = -211.8087,
    c2 = -1093451,
    w_ref = 439110,
    z = -1,
    a_0 = 3.72,
    MM = 0.18925) "GdO2-";

  //943
  constant DataRecord GdHCO3_p2(
    R = Modelica.Constants.R/GdHCO3_p2.MM,
    G_ref = -1260221,
    H_ref = -1366913,
    S_ref = -40.585,
    a1 = 2.5213e-06,
    a2 = -2637.72,
    a3 = 0.00034373,
    a4 = -105366,
    c1 = 163.7601,
    c2 = 191129,
    w_ref = 504090,
    z = 2,
    a_0 = 3.72,
    MM = 0.21827) "GdHCO3+2";

  //944
  constant DataRecord GdCl_p2(
    R = Modelica.Constants.R/GdCl_p2.MM,
    G_ref = -796634,
    H_ref = -839310,
    S_ref = -93.303,
    a1 = -1.1004e-06,
    a2 = -3521.67,
    a3 = 0.00037834,
    a4 = -101709,
    c1 = 75.8275,
    c2 = -140733,
    w_ref = 586010,
    z = 2,
    a_0 = 3.72,
    MM = 0.1927) "GdCl+2";

  //945
  constant DataRecord GdCl2_p1(
    R = Modelica.Constants.R/GdCl2_p1.MM,
    G_ref = -925919,
    H_ref = -999976,
    S_ref = -22.594,
    a1 = 1.1921e-05,
    a2 = -346.1,
    a3 = 0.0002544,
    a4 = -114838,
    c1 = 39.8568,
    c2 = -162582,
    w_ref = 263800,
    z = 1,
    a_0 = 3.72,
    MM = 0.22815) "GdCl2+";

  //946
  constant DataRecord GdCl3_aq(
    R = Modelica.Constants.R/GdCl3_aq.MM,
    G_ref = -1054786,
    H_ref = -1172231,
    S_ref = 9.205,
    a1 = 2.6709e-05,
    a2 = 3264.69,
    a3 = 0.0001125,
    a4 = -129767,
    c1 = -82.2428,
    c2 = -498473,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15725) "GdCl3,aq";

  //947
  constant DataRecord GdCl4_n1(
    R = Modelica.Constants.R/GdCl4_n1.MM,
    G_ref = -1183654,
    H_ref = -1356871,
    S_ref = -0.837,
    a1 = 4.6575e-05,
    a2 = 8116.75,
    a3 = -7.8496e-05,
    a4 = -149825,
    c1 = -205.2114,
    c2 = -1148416,
    w_ref = 682410,
    z = -1,
    a_0 = 3.72,
    MM = 0.29905) "GdCl4-";

  //948
  constant DataRecord GdH2PO4_p2(
    R = Modelica.Constants.R/GdH2PO4_p2.MM,
    G_ref = -1806651,
    H_ref = -1994094,
    S_ref = -109.202,
    a1 = 6.7145e-06,
    a2 = -1616.36,
    a3 = 0.00030412,
    a4 = -109587,
    c1 = 180.1417,
    c2 = 214225,
    w_ref = 609780,
    z = 2,
    a_0 = 3.72,
    MM = 0.25424) "GdH2PO4+2";

  //949
  constant DataRecord GdNO3_p2(
    R = Modelica.Constants.R/GdNO3_p2.MM,
    G_ref = -776969,
    H_ref = -919643,
    S_ref = -135.98,
    a1 = 5.525e-06,
    a2 = -1905.18,
    a3 = 0.00031515,
    a4 = -108395,
    c1 = 140.8372,
    c2 = 65534,
    w_ref = 647470,
    z = 2,
    a_0 = 3.72,
    MM = 0.21926) "GdNO3+2";

  //950
  constant DataRecord GdF_p2(
    R = Modelica.Constants.R/GdF_p2.MM,
    G_ref = -971525,
    H_ref = -1001231,
    S_ref = -60.25,
    a1 = -1.2216e-05,
    a2 = -6235.88,
    a3 = 0.00048507,
    a4 = -90492,
    c1 = 81.1432,
    c2 = -105776,
    w_ref = 534550,
    z = 2,
    a_0 = 3.72,
    MM = 0.17625) "GdF+2";

  //951
  constant DataRecord GdF2_p1(
    R = Modelica.Constants.R/GdF2_p1.MM,
    G_ref = -1272354,
    H_ref = -1346411,
    S_ref = -42.677,
    a1 = -1.1098e-05,
    a2 = -5964.58,
    a3 = 0.0004747,
    a4 = -91613,
    c1 = 66.9239,
    c2 = -79103,
    w_ref = 296900,
    z = 1,
    a_0 = 3.72,
    MM = 0.19525) "GdF2+";

  //952
  constant DataRecord GdF3_aq(
    R = Modelica.Constants.R/GdF3_aq.MM,
    G_ref = -1568582,
    H_ref = -1704562,
    S_ref = -83.262,
    a1 = -1.0005e-05,
    a2 = -5697.52,
    a3 = 0.00046421,
    a4 = -92717,
    c1 = -40.3635,
    c2 = -352916,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15725) "GdF3,aq";

  //953
  constant DataRecord GdF4_n1(
    R = Modelica.Constants.R/GdF4_n1.MM,
    G_ref = -1861880,
    H_ref = -2080703,
    S_ref = -194.138,
    a1 = -4.2928e-06,
    a2 = -4301.36,
    a3 = 0.00040901,
    a4 = -98487,
    c1 = -114.8587,
    c2 = -927212,
    w_ref = 972320,
    z = -1,
    a_0 = 3.72,
    MM = 0.23325) "GdF4-";

  //954
  constant DataRecord GdSO4_p1(
    R = Modelica.Constants.R/GdSO4_p1.MM,
    G_ref = -1386996,
    H_ref = -1576531,
    S_ref = -50.208,
    a1 = -6.1823e-06,
    a2 = -1744.94,
    a3 = 0.00030887,
    a4 = -109056,
    c1 = -87.5962,
    c2 = -618730,
    w_ref = 304890,
    z = 1,
    a_0 = 3.72,
    MM = 0.25331) "GdSO4+";

  //955
  constant DataRecord TbCO3_p1(
    R = Modelica.Constants.R/TbCO3_p1.MM,
    G_ref = -1240556,
    H_ref = -1394653,
    S_ref = -195.393,
    a1 = 6.358e-06,
    a2 = -1700.5,
    a3 = 0.00030679,
    a4 = -109240,
    c1 = -85.0833,
    c2 = -619583,
    w_ref = 334800,
    z = 1,
    a_0 = 3.72,
    MM = 0.21893) "TbCO3+";

  //956
  constant DataRecord TbOH_p2(
    R = Modelica.Constants.R/TbOH_p2.MM,
    G_ref = -859812,
    H_ref = -904372,
    S_ref = -38.911,
    a1 = 1.1024e-05,
    a2 = -562.79,
    a3 = 0.00026242,
    a4 = -113943,
    c1 = 6.7124,
    c2 = -353669,
    w_ref = 500780,
    z = 2,
    a_0 = 3.72,
    MM = 0.17593) "TbOH+2";

  //957
  constant DataRecord TbO_p1(
    R = Modelica.Constants.R/TbO_p1.MM,
    G_ref = -812114,
    H_ref = -874456,
    S_ref = 26.359,
    a1 = 1.1407e-05,
    a2 = -469.86,
    a3 = 0.00025889,
    a4 = -114328,
    c1 = -123.3841,
    c2 = -706514,
    w_ref = 190580,
    z = 1,
    a_0 = 3.72,
    MM = 0.17492) "TbO+";

  //958
  constant DataRecord TbO2H_aq(
    R = Modelica.Constants.R/TbO2H_aq.MM,
    G_ref = -998721,
    H_ref = -1050184,
    S_ref = 171.126,
    a1 = 1.9881e-05,
    a2 = 1599.46,
    a3 = 0.00017751,
    a4 = -122880,
    c1 = -262.9845,
    c2 = -1126688,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19092) "TbO2H,aq";

  //959
  constant DataRecord TbO2_n1(
    R = Modelica.Constants.R/TbO2_n1.MM,
    G_ref = -946421,
    H_ref = -1010854,
    S_ref = 127.612,
    a1 = 2.0086e-05,
    a2 = 1647.66,
    a3 = 0.00017603,
    a4 = -123081,
    c1 = -168.0253,
    c2 = -957086,
    w_ref = 488520,
    z = -1,
    a_0 = 3.72,
    MM = 0.19092) "TbO2-";

  //960
  constant DataRecord TbHCO3_p2(
    R = Modelica.Constants.R/TbHCO3_p2.MM,
    G_ref = -1263986,
    H_ref = -1402895,
    S_ref = -143.511,
    a1 = 3.1782e-06,
    a2 = -2478.27,
    a3 = 0.00033765,
    a4 = -106023,
    c1 = 177.6351,
    c2 = 189213,
    w_ref = 660650,
    z = 2,
    a_0 = 3.72,
    MM = 0.21994) "TbHCO3+2";

  //961
  constant DataRecord TbCl_p2(
    R = Modelica.Constants.R/TbCl_p2.MM,
    G_ref = -799981,
    H_ref = -851444,
    S_ref = -117.57,
    a1 = -8.51e-07,
    a2 = -3461.67,
    a3 = 0.00037618,
    a4 = -101960,
    c1 = 78.5943,
    c2 = -142649,
    w_ref = 622040,
    z = 2,
    a_0 = 3.72,
    MM = 0.19437) "TbCl+2";

  //962
  constant DataRecord TbCl2_p1(
    R = Modelica.Constants.R/TbCl2_p1.MM,
    G_ref = -929685,
    H_ref = -1014202,
    S_ref = -52.718,
    a1 = 1.2216e-05,
    a2 = -271.88,
    a3 = 0.00025101,
    a4 = -115144,
    c1 = 42.9404,
    c2 = -166322,
    w_ref = 308950,
    z = 1,
    a_0 = 3.72,
    MM = 0.22982) "TbCl2+";

  //963
  constant DataRecord TbCl3_aq(
    R = Modelica.Constants.R/TbCl3_aq.MM,
    G_ref = -1058552,
    H_ref = -1189511,
    S_ref = -30.125,
    a1 = 2.6867e-05,
    a2 = 3304.19,
    a3 = 0.00011079,
    a4 = -129930,
    c1 = -84.062,
    c2 = -504800,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15892) "TbCl3,aq";

  //964
  constant DataRecord TbCl4_n1(
    R = Modelica.Constants.R/TbCl4_n1.MM,
    G_ref = -1187419,
    H_ref = -1378210,
    S_ref = -55.229,
    a1 = 4.7027e-05,
    a2 = 8228.13,
    a3 = -8.3098e-05,
    a4 = -150285,
    c1 = -200.4851,
    c2 = -1158081,
    w_ref = 763870,
    z = -1,
    a_0 = 3.72,
    MM = 0.30072) "TbCl4-";

  //965
  constant DataRecord TbH2PO4_p2(
    R = Modelica.Constants.R/TbH2PO4_p2.MM,
    G_ref = -1809998,
    H_ref = -2007902,
    S_ref = -138.49,
    a1 = 6.9844e-06,
    a2 = -1551.3,
    a3 = 0.00030173,
    a4 = -109855,
    c1 = 183.4671,
    c2 = 212309,
    w_ref = 651830,
    z = 2,
    a_0 = 3.72,
    MM = 0.25591) "TbH2PO4+2";

  //966
  constant DataRecord TbNO3_p2(
    R = Modelica.Constants.R/TbNO3_p2.MM,
    G_ref = -781153,
    H_ref = -936379,
    S_ref = -174.054,
    a1 = 5.8538e-06,
    a2 = -1824.27,
    a3 = 0.00031181,
    a4 = -108730,
    c1 = 145.766,
    c2 = 63618,
    w_ref = 706970,
    z = 2,
    a_0 = 3.72,
    MM = 0.22093) "TbNO3+2";

  //967
  constant DataRecord TbF_p2(
    R = Modelica.Constants.R/TbF_p2.MM,
    G_ref = -975709,
    H_ref = -1010854,
    S_ref = -73.638,
    a1 = -1.2016e-05,
    a2 = -6186.42,
    a3 = 0.00048299,
    a4 = -90697,
    c1 = 82.5645,
    c2 = -107692,
    w_ref = 555970,
    z = 2,
    a_0 = 3.72,
    MM = 0.17792) "TbF+2";

  //968
  constant DataRecord TbF2_p1(
    R = Modelica.Constants.R/TbF2_p1.MM,
    G_ref = -1277375,
    H_ref = -1356871,
    S_ref = -55.647,
    a1 = -1.0902e-05,
    a2 = -5914.29,
    a3 = 0.00047226,
    a4 = -91818,
    c1 = 67.339,
    c2 = -82843,
    w_ref = 313090,
    z = 1,
    a_0 = 3.72,
    MM = 0.19692) "TbF2+";

  //969
  constant DataRecord TbF3_aq(
    R = Modelica.Constants.R/TbF3_aq.MM,
    G_ref = -1573602,
    H_ref = -1716277,
    S_ref = -100.416,
    a1 = -9.8466e-06,
    a2 = -5658.15,
    a3 = 0.00046248,
    a4 = -92876,
    c1 = -42.1831,
    c2 = -359238,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15892) "TbF3,aq";

  //970
  constant DataRecord TbF4_n1(
    R = Modelica.Constants.R/TbF4_n1.MM,
    G_ref = -1867319,
    H_ref = -2095766,
    S_ref = -222.589,
    a1 = -3.9731e-06,
    a2 = -4223.25,
    a3 = 0.00040595,
    a4 = -98809,
    c1 = -113.7291,
    c2 = -936877,
    w_ref = 1014790,
    z = -1,
    a_0 = 3.72,
    MM = 0.23492) "TbF4-";

  //971
  constant DataRecord TbSO4_p1(
    R = Modelica.Constants.R/TbSO4_p1.MM,
    G_ref = -1432602,
    H_ref = -1588246,
    S_ref = -71.546,
    a1 = 6.3563e-06,
    a2 = -1702.51,
    a3 = 0.00030721,
    a4 = -109232,
    c1 = -84.6641,
    c2 = -619583,
    w_ref = 339360,
    z = 1,
    a_0 = 3.72,
    MM = 0.25498) "TbSO4+";

  //972
  constant DataRecord DyCO3_p1(
    R = Modelica.Constants.R/DyCO3_p1.MM,
    G_ref = -1237627,
    H_ref = -1393941,
    S_ref = -201.25,
    a1 = -3.9079e-06,
    a2 = -4206.43,
    a3 = 0.0004051,
    a4 = -98880,
    c1 = -68.3293,
    c2 = -626399,
    w_ref = 537980,
    z = 1,
    a_0 = 3.72,
    MM = 0.22251) "DyCO3+";

  //973
  constant DataRecord DyHCO3_p2(
    R = Modelica.Constants.R/DyHCO3_p2.MM,
    G_ref = -1260639,
    H_ref = -1379465,
    S_ref = -76.149,
    a1 = 2.5175e-06,
    a2 = -2639.35,
    a3 = 0.00034395,
    a4 = -105357,
    c1 = 192.0502,
    c2 = 271671,
    w_ref = 559610,
    z = 2,
    a_0 = 3.72,
    MM = 0.22352) "DyHCO3+2";

  //974
  constant DataRecord DyCl_p2(
    R = Modelica.Constants.R/DyCl_p2.MM,
    G_ref = -796634,
    H_ref = -850189,
    S_ref = -123.428,
    a1 = -1.1418e-06,
    a2 = -3531.71,
    a3 = 0.00037872,
    a4 = -101667,
    c1 = 103.0875,
    c2 = -60191,
    w_ref = 630400,
    z = 2,
    a_0 = 3.72,
    MM = 0.19795) "DyCl+2";

  //975
  constant DataRecord DyCl2_p1(
    R = Modelica.Constants.R/DyCl2_p1.MM,
    G_ref = -926338,
    H_ref = -1013365,
    S_ref = -60.25,
    a1 = 1.1903e-05,
    a2 = -349.24,
    a3 = 0.00025426,
    a4 = -114826,
    c1 = 90.3916,
    c2 = -5439,
    w_ref = 321580,
    z = 1,
    a_0 = 3.72,
    MM = 0.2334) "DyCl2+";

  //976
  constant DataRecord DyCl3_aq(
    R = Modelica.Constants.R/DyCl3_aq.MM,
    G_ref = -1055205,
    H_ref = -1189093,
    S_ref = -40.166,
    a1 = 2.6471e-05,
    a2 = 3206.87,
    a3 = 0.00011474,
    a4 = -129528,
    c1 = -5.8262,
    c2 = -232869,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1625) "DyCl3,aq";

  //977
  constant DataRecord DyCl4_n1(
    R = Modelica.Constants.R/DyCl4_n1.MM,
    G_ref = -1184072,
    H_ref = -1379046,
    S_ref = -68.618,
    a1 = 4.6659e-05,
    a2 = 8136.5,
    a3 = -7.9115e-05,
    a4 = -149904,
    c1 = -78.9161,
    c2 = -742493,
    w_ref = 785590,
    z = -1,
    a_0 = 3.72,
    MM = 0.3043) "DyCl4-";

  //978
  constant DataRecord DyH2PO4_p2(
    R = Modelica.Constants.R/DyH2PO4_p2.MM,
    G_ref = -1806651,
    H_ref = -2007065,
    S_ref = -145.603,
    a1 = 6.7103e-06,
    a2 = -1617.87,
    a3 = 0.00030427,
    a4 = -109579,
    c1 = 208.4147,
    c2 = 294767,
    w_ref = 665130,
    z = 2,
    a_0 = 3.72,
    MM = 0.25949) "DyH2PO4+2";

  //979
  constant DataRecord DyNO3_p2(
    R = Modelica.Constants.R/DyNO3_p2.MM,
    G_ref = -775714,
    H_ref = -933869,
    S_ref = -183.259,
    a1 = 5.5844e-06,
    a2 = -1889.45,
    a3 = 0.00031424,
    a4 = -108458,
    c1 = 170.8398,
    c2 = 146076,
    w_ref = 721610,
    z = 2,
    a_0 = 3.72,
    MM = 0.22451) "DyNO3+2";

  //980
  constant DataRecord DyF_p2(
    R = Modelica.Constants.R/DyF_p2.MM,
    G_ref = -972362,
    H_ref = -1008762,
    S_ref = -76.986,
    a1 = -1.2322e-05,
    a2 = -6262.74,
    a3 = 0.0004863,
    a4 = -90379,
    c1 = 106.6259,
    c2 = -25234,
    w_ref = 559610,
    z = 2,
    a_0 = 3.72,
    MM = 0.1815) "DyF+2";

  //981
  constant DataRecord DyF2_p1(
    R = Modelica.Constants.R/DyF2_p1.MM,
    G_ref = -1274028,
    H_ref = -1354779,
    S_ref = -58.994,
    a1 = -1.1229e-05,
    a2 = -5995.76,
    a3 = 0.00047583,
    a4 = -91483,
    c1 = 114.4098,
    c2 = 78036,
    w_ref = 321580,
    z = 1,
    a_0 = 3.72,
    MM = 0.2005) "DyF2+";

  //982
  constant DataRecord DyF3_aq(
    R = Modelica.Constants.R/DyF3_aq.MM,
    G_ref = -1570674,
    H_ref = -1714603,
    S_ref = -104.6,
    a1 = -1.0243e-05,
    a2 = -5755.38,
    a3 = 0.00046645,
    a4 = -92475,
    c1 = 36.0527,
    c2 = -87312,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1625) "DyF3,aq";

  //983
  constant DataRecord DyF4_n1(
    R = Modelica.Constants.R/DyF4_n1.MM,
    G_ref = -1864390,
    H_ref = -2095347,
    S_ref = -229.702,
    a1 = -4.3522e-06,
    a2 = -4317.13,
    a3 = 0.00040992,
    a4 = -98420,
    c1 = 7.5304,
    c2 = -521289,
    w_ref = 1033110,
    z = -1,
    a_0 = 3.72,
    MM = 0.2385) "DyF4-";

  //984
  constant DataRecord DyOH_p2(
    R = Modelica.Constants.R/DyOH_p2.MM,
    G_ref = -856465,
    H_ref = -903493,
    S_ref = -45.606,
    a1 = 1.0943e-05,
    a2 = -583.29,
    a3 = 0.00026338,
    a4 = -113859,
    c1 = 12.0369,
    c2 = -338331,
    w_ref = 510700,
    z = 2,
    a_0 = 3.72,
    MM = 0.17951) "DyOH+2";

  //985
  constant DataRecord DyO_p1(
    R = Modelica.Constants.R/DyO_p1.MM,
    G_ref = -809186,
    H_ref = -836884,
    S_ref = 19.246,
    a1 = 1.127e-05,
    a2 = -504.42,
    a3 = 0.00026048,
    a4 = -114186,
    c1 = -117.7842,
    c2 = -690322,
    w_ref = 200790,
    z = 1,
    a_0 = 3.72,
    MM = 0.1785) "DyO+";

  //986
  constant DataRecord DyO2H_aq(
    R = Modelica.Constants.R/DyO2H_aq.MM,
    G_ref = -996629,
    H_ref = -1050602,
    S_ref = 163.594,
    a1 = 1.9652e-05,
    a2 = 1542.01,
    a3 = 0.00018013,
    a4 = -122646,
    c1 = -252.6856,
    c2 = -1090894,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1945) "DyO2H,aq";

  //987
  constant DataRecord DyO2_n1(
    R = Modelica.Constants.R/DyO2_n1.MM,
    G_ref = -947258,
    H_ref = -1014536,
    S_ref = 119.662,
    a1 = 1.9838e-05,
    a2 = 1587.24,
    a3 = 0.00017839,
    a4 = -122830,
    c1 = -157.1632,
    c2 = -922995,
    w_ref = 499950,
    z = -1,
    a_0 = 3.72,
    MM = 0.1945) "DyO2-";

  //988
  constant DataRecord DySO4_p1(
    R = Modelica.Constants.R/DySO4_p1.MM,
    G_ref = -1429254,
    H_ref = -1585736,
    S_ref = -74.894,
    a1 = -6.2576e-06,
    a2 = -1726.65,
    a3 = 0.00030815,
    a4 = -109131,
    c1 = -73.6936,
    c2 = -582936,
    w_ref = 344010,
    z = 1,
    a_0 = 3.72,
    MM = 0.25856) "DySO4+";

  //989
  constant DataRecord HoCO3_p1(
    R = Modelica.Constants.R/HoCO3_p1.MM,
    G_ref = -1249342,
    H_ref = -1404318,
    S_ref = -196.648,
    a1 = -4.1648e-06,
    a2 = -4269.65,
    a3 = 0.00040766,
    a4 = -98617,
    c1 = -73.0175,
    c2 = -640035,
    w_ref = 529650,
    z = 1,
    a_0 = 3.72,
    MM = 0.22494) "HoCO3+";

  //990
  constant DataRecord HoHCO3_p2(
    R = Modelica.Constants.R/HoHCO3_p2.MM,
    G_ref = -1271936,
    H_ref = -1389506,
    S_ref = -71.128,
    a1 = 1.9066e-06,
    a2 = -2787.38,
    a3 = 0.00034952,
    a4 = -104746,
    c1 = 182.2182,
    c2 = 240990,
    w_ref = 548730,
    z = 2,
    a_0 = 3.72,
    MM = 0.22595) "HoHCO3+2";

  //991
  constant DataRecord HoCl_p2(
    R = Modelica.Constants.R/HoCl_p2.MM,
    G_ref = -807930,
    H_ref = -860230,
    S_ref = -118.826,
    a1 = -1.7439e-06,
    a2 = -3679.41,
    a3 = 0.00038465,
    a4 = -101060,
    c1 = 93.4915,
    c2 = -90876,
    w_ref = 622040,
    z = 2,
    a_0 = 3.72,
    MM = 0.20038) "HoCl+2";

  //992
  constant DataRecord HoCl2_p1(
    R = Modelica.Constants.R/HoCl2_p1.MM,
    G_ref = -937634,
    H_ref = -1023406,
    S_ref = -54.81,
    a1 = 1.1235e-05,
    a2 = -513.46,
    a3 = 0.00026096,
    a4 = -114148,
    c1 = 72.3861,
    c2 = -65300,
    w_ref = 313090,
    z = 1,
    a_0 = 3.72,
    MM = 0.23583) "HoCl2+";

  //993
  constant DataRecord HoCl3_aq(
    R = Modelica.Constants.R/HoCl3_aq.MM,
    G_ref = -1066502,
    H_ref = -1198298,
    S_ref = -32.635,
    a1 = 2.5759e-05,
    a2 = 3033.23,
    a3 = 0.00012146,
    a4 = -128809,
    c1 = -34.9364,
    c2 = -334055,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16493) "HoCl3,aq";

  //994
  constant DataRecord HoCl4_n1(
    R = Modelica.Constants.R/HoCl4_n1.MM,
    G_ref = -1195369,
    H_ref = -1387833,
    S_ref = -58.576,
    a1 = 4.582e-05,
    a2 = 7933.32,
    a3 = -7.1467e-05,
    a4 = -149068,
    c1 = -124.6259,
    c2 = -897129,
    w_ref = 772370,
    z = -1,
    a_0 = 3.72,
    MM = 0.30673) "HoCl4-";

  //995
  constant DataRecord HoH2PO4_p2(
    R = Modelica.Constants.R/HoH2PO4_p2.MM,
    G_ref = -1818366,
    H_ref = -2017106,
    S_ref = -140.164,
    a1 = 6.1061e-06,
    a2 = -1762.43,
    a3 = 0.00030932,
    a4 = -108985,
    c1 = 198.768,
    c2 = 264086,
    w_ref = 656220,
    z = 2,
    a_0 = 3.72,
    MM = 0.26192) "HoH2PO4+2";

  //996
  constant DataRecord HoNO3_p2(
    R = Modelica.Constants.R/HoNO3_p2.MM,
    G_ref = -787429,
    H_ref = -943910,
    S_ref = -176.146,
    a1 = 4.9769e-06,
    a2 = -2038.32,
    a3 = 0.00032023,
    a4 = -107843,
    c1 = 161.1091,
    c2 = 115391,
    w_ref = 711820,
    z = 2,
    a_0 = 3.72,
    MM = 0.22694) "HoNO3+2";

  //997
  constant DataRecord HoF_p2(
    R = Modelica.Constants.R/HoF_p2.MM,
    G_ref = -984077,
    H_ref = -1020059,
    S_ref = -74.475,
    a1 = -1.2909e-05,
    a2 = -6404.2,
    a3 = 0.00049146,
    a4 = -89793,
    c1 = 97.4608,
    c2 = -55915,
    w_ref = 555970,
    z = 2,
    a_0 = 3.72,
    MM = 0.18393) "HoF+2";

  //998
  constant DataRecord HoF2_p1(
    R = Modelica.Constants.R/HoF2_p1.MM,
    G_ref = -1286162,
    H_ref = -1366076,
    S_ref = -56.484,
    a1 = -1.1883e-05,
    a2 = -6153.41,
    a3 = 0.00048157,
    a4 = -90830,
    c1 = 96.7927,
    c2 = 18171,
    w_ref = 317310,
    z = 1,
    a_0 = 3.72,
    MM = 0.20293) "HoF2+";

  //999
  constant DataRecord HoF3_aq(
    R = Modelica.Constants.R/HoF3_aq.MM,
    G_ref = -1582807,
    H_ref = -1725900,
    S_ref = -101.253,
    a1 = -1.0955e-05,
    a2 = -5928.98,
    a3 = 0.00047317,
    a4 = -91759,
    c1 = 6.9417,
    c2 = -188493,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16493) "HoF3,aq";

  //1000
  constant DataRecord HoF4_n1(
    R = Modelica.Constants.R/HoF4_n1.MM,
    G_ref = -1876524,
    H_ref = -2106226,
    S_ref = -224.681,
    a1 = -5.1777e-06,
    a2 = -4518.59,
    a3 = 0.00041781,
    a4 = -97588,
    c1 = -37.8158,
    c2 = -675925,
    w_ref = 1023870,
    z = -1,
    a_0 = 3.72,
    MM = 0.24093) "HoF4-";

  //1001
  constant DataRecord HoOH_p2(
    R = Modelica.Constants.R/HoOH_p2.MM,
    G_ref = -868180,
    H_ref = -913786,
    S_ref = -40.585,
    a1 = 1.1035e-05,
    a2 = -560.4,
    a3 = 0.0002624,
    a4 = -113951,
    c1 = 7.9948,
    c2 = -350259,
    w_ref = 504090,
    z = 2,
    a_0 = 3.72,
    MM = 0.18194) "HoOH+2";

  //1002
  constant DataRecord HoO_p1(
    R = Modelica.Constants.R/HoO_p1.MM,
    G_ref = -820901,
    H_ref = -847009,
    S_ref = 24.686,
    a1 = 1.1358e-05,
    a2 = -483.17,
    a3 = 0.00025973,
    a4 = -114273,
    c1 = -121.9268,
    c2 = -702255,
    w_ref = 193090,
    z = 1,
    a_0 = 3.72,
    MM = 0.18093) "HoO+";

  //1003
  constant DataRecord HoO2H_aq(
    R = Modelica.Constants.R/HoO2H_aq.MM,
    G_ref = -1009599,
    H_ref = -1061899,
    S_ref = 169.452,
    a1 = 1.9824e-05,
    a2 = 1585.74,
    a3 = 0.000178,
    a4 = -122826,
    c1 = -260.5327,
    c2 = -1118166,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19693) "HoO2H,aq";

  //1004
  constant DataRecord HoO2_n1(
    R = Modelica.Constants.R/HoO2_n1.MM,
    G_ref = -958554,
    H_ref = -1024160,
    S_ref = 125.52,
    a1 = 2.0039e-05,
    a2 = 1636.45,
    a3 = 0.0001764,
    a4 = -123035,
    c1 = -165.2948,
    c2 = -948563,
    w_ref = 491540,
    z = -1,
    a_0 = 3.72,
    MM = 0.19693) "HoO2-";

  //1005
  constant DataRecord HoSO4_p1(
    R = Modelica.Constants.R/HoSO4_p1.MM,
    G_ref = -1440133,
    H_ref = -1596196,
    S_ref = -71.128,
    a1 = -6.0128e-06,
    a2 = -1786.4,
    a3 = 0.0003105,
    a4 = -108884,
    c1 = -78.0433,
    c2 = -596571,
    w_ref = 339360,
    z = 1,
    a_0 = 3.72,
    MM = 0.26099) "HoSO4+";

  //1006
  constant DataRecord ErCO3_p1(
    R = Modelica.Constants.R/ErCO3_p1.MM,
    G_ref = -1243485,
    H_ref = -1404025,
    S_ref = -217.15,
    a1 = -4.4095e-06,
    a2 = -4330.86,
    a3 = 0.00041039,
    a4 = -98366,
    c1 = -72.3075,
    c2 = -648558,
    w_ref = 564000,
    z = 1,
    a_0 = 3.72,
    MM = 0.22727) "ErCO3+";

  //1007
  constant DataRecord ErHCO3_p2(
    R = Modelica.Constants.R/ErHCO3_p2.MM,
    G_ref = -1266078,
    H_ref = -1389925,
    S_ref = -94.558,
    a1 = 1.1397e-06,
    a2 = -2975.49,
    a3 = 0.00035708,
    a4 = -103968,
    c1 = 180.1367,
    c2 = 221815,
    w_ref = 586010,
    z = 2,
    a_0 = 3.72,
    MM = 0.22828) "ErHCO3+2";

  //1008
  constant DataRecord ErCl_p2(
    R = Modelica.Constants.R/ErCl_p2.MM,
    G_ref = -802073,
    H_ref = -859394,
    S_ref = -138.909,
    a1 = -2.5363e-06,
    a2 = -3873.88,
    a3 = 0.0003925,
    a4 = -100253,
    c1 = 90.7196,
    c2 = -110052,
    w_ref = 651830,
    z = 2,
    a_0 = 3.72,
    MM = 0.20271) "ErCl+2";

  //1009
  constant DataRecord ErCl2_p1(
    R = Modelica.Constants.R/ErCl2_p1.MM,
    G_ref = -931358,
    H_ref = -1023825,
    S_ref = -79.496,
    a1 = 1.0376e-05,
    a2 = -721.95,
    a3 = 0.00026891,
    a4 = -113286,
    c1 = 65.3428,
    c2 = -102717,
    w_ref = 353510,
    z = 1,
    a_0 = 3.72,
    MM = 0.23816) "ErCl2+";

  //1010
  constant DataRecord ErCl3_aq(
    R = Modelica.Constants.R/ErCl3_aq.MM,
    G_ref = -1060226,
    H_ref = -1201226,
    S_ref = -65.27,
    a1 = 2.465e-05,
    a2 = 2762.4,
    a3 = 0.00013213,
    a4 = -127687,
    c1 = -53.1314,
    c2 = -397292,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16726) "ErCl3,aq";

  //1011
  constant DataRecord ErCl4_n1(
    R = Modelica.Constants.R/ErCl4_n1.MM,
    G_ref = -1189093,
    H_ref = -1394109,
    S_ref = -103.345,
    a1 = 4.4814e-05,
    a2 = 7687.85,
    a3 = -6.1877e-05,
    a4 = -148051,
    c1 = -146.189,
    c2 = -993779,
    w_ref = 840150,
    z = -1,
    a_0 = 3.72,
    MM = 0.30906) "ErCl4-";

  //1012
  constant DataRecord ErH2PO4_p2(
    R = Modelica.Constants.R/ErH2PO4_p2.MM,
    G_ref = -1812090,
    H_ref = -2017525,
    S_ref = -164.431,
    a1 = 5.3363e-06,
    a2 = -1949.74,
    a3 = 0.00031654,
    a4 = -108211,
    c1 = 196.6091,
    c2 = 244906,
    w_ref = 692700,
    z = 2,
    a_0 = 3.72,
    MM = 0.26425) "ErH2PO4+2";

  //1013
  constant DataRecord ErNO3_p2(
    R = Modelica.Constants.R/ErNO3_p2.MM,
    G_ref = -780734,
    H_ref = -945584,
    S_ref = -207.945,
    a1 = 4.2376e-06,
    a2 = -2217.98,
    a3 = 0.00032711,
    a4 = -107102,
    c1 = 159.7832,
    c2 = 96215,
    w_ref = 757300,
    z = 2,
    a_0 = 3.72,
    MM = 0.22927) "ErNO3+2";

  //1014
  constant DataRecord ErF_p2(
    R = Modelica.Constants.R/ErF_p2.MM,
    G_ref = -977801,
    H_ref = -1016294,
    S_ref = -85.354,
    a1 = -1.3752e-05,
    a2 = -6611.47,
    a3 = 0.00049993,
    a4 = -88939,
    c1 = 93.3082,
    c2 = -75094,
    w_ref = 570780,
    z = 2,
    a_0 = 3.72,
    MM = 0.18626) "ErF+2";

  //1015
  constant DataRecord ErF2_p1(
    R = Modelica.Constants.R/ErF2_p1.MM,
    G_ref = -1279886,
    H_ref = -1362729,
    S_ref = -67.362,
    a1 = -1.2818e-05,
    a2 = -6383.36,
    a3 = 0.00049092,
    a4 = -89881,
    c1 = 87.6397,
    c2 = -19238,
    w_ref = 334800,
    z = 1,
    a_0 = 3.72,
    MM = 0.20526) "ErF2+";

  //1016
  constant DataRecord ErF3_aq(
    R = Modelica.Constants.R/ErF3_aq.MM,
    G_ref = -1576531,
    H_ref = -1723390,
    S_ref = -115.897,
    a1 = -1.2064e-05,
    a2 = -6199.81,
    a3 = 0.00048384,
    a4 = -90638,
    c1 = -11.252,
    c2 = -251735,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16726) "ErF3,aq";

  //1017
  constant DataRecord ErF4_n1(
    R = Modelica.Constants.R/ErF4_n1.MM,
    G_ref = -1870666,
    H_ref = -2106644,
    S_ref = -248.111,
    a1 = -6.2827e-06,
    a2 = -4787.79,
    a3 = 0.00042823,
    a4 = -96475,
    c1 = -62.0768,
    c2 = -772576,
    w_ref = 1062320,
    z = -1,
    a_0 = 3.72,
    MM = 0.24326) "ErF4-";

  //1018
  constant DataRecord ErOH_p2(
    R = Modelica.Constants.R/ErOH_p2.MM,
    G_ref = -861904,
    H_ref = -913032,
    S_ref = -61.086,
    a1 = 1.068e-05,
    a2 = -647.6,
    a3 = 0.00026592,
    a4 = -113591,
    c1 = 24.5337,
    c2 = -302532,
    w_ref = 534550,
    z = 2,
    a_0 = 3.72,
    MM = 0.18427) "ErOH+2";

  //1019
  constant DataRecord ErO_p1(
    R = Modelica.Constants.R/ErO_p1.MM,
    G_ref = -815043,
    H_ref = -847218,
    S_ref = 2.51,
    a1 = 1.101e-05,
    a2 = -566.97,
    a3 = 0.00026277,
    a4 = -113926,
    c1 = -104.2101,
    c2 = -651118,
    w_ref = 225680,
    z = 1,
    a_0 = 3.72,
    MM = 0.18326) "ErO+";

  //1020
  constant DataRecord ErO2H_aq(
    R = Modelica.Constants.R/ErO2H_aq.MM,
    G_ref = -1004578,
    H_ref = -1063991,
    S_ref = 145.603,
    a1 = 1.9022e-05,
    a2 = 1388.5,
    a3 = 0.00018609,
    a4 = -122010,
    c1 = -228.9004,
    c2 = -1008223,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19926) "ErO2H,aq";

  //1021
  constant DataRecord ErO2_n1(
    R = Modelica.Constants.R/ErO2_n1.MM,
    G_ref = -957299,
    H_ref = -1029808,
    S_ref = 100.416,
    a1 = 1.931e-05,
    a2 = 1459.25,
    a3 = 0.00018323,
    a4 = -122303,
    c1 = -131.3399,
    c2 = -842879,
    w_ref = 530070,
    z = -1,
    a_0 = 3.72,
    MM = 0.19926) "ErO2-";

  //1022
  constant DataRecord ErSO4_p1(
    R = Modelica.Constants.R/ErSO4_p1.MM,
    G_ref = -1433857,
    H_ref = -1594305,
    S_ref = -88.282,
    a1 = 5.7501e-06,
    a2 = -1850.5,
    a3 = 0.00031302,
    a4 = -108621,
    c1 = -78.2906,
    c2 = -605094,
    w_ref = 363300,
    z = 1,
    a_0 = 3.72,
    MM = 0.26332) "ErSO4+";

  //1023
  constant DataRecord TmCO3_p1(
    R = Modelica.Constants.R/TmCO3_p1.MM,
    G_ref = -1243903,
    H_ref = -1308337,
    S_ref = -215.894,
    a1 = -4.4396e-06,
    a2 = -4338.6,
    a3 = 0.00041079,
    a4 = -98332,
    c1 = -73.1246,
    c2 = -648558,
    w_ref = 555130,
    z = 1,
    a_0 = 3.72,
    MM = 0.22894) "TmCO3+";

  //1024
  constant DataRecord TmHCO3_p2(
    R = Modelica.Constants.R/TmHCO3_p2.MM,
    G_ref = -1266078,
    H_ref = -1389925,
    S_ref = -93.722,
    a1 = 1.1397e-06,
    a2 = -2975.49,
    a3 = 0.00035708,
    a4 = -103968,
    c1 = 180.1367,
    c2 = 221815,
    w_ref = 586010,
    z = 2,
    a_0 = 3.72,
    MM = 0.22995) "TmHCO3+2";

  //1025
  constant DataRecord TmCl_p2(
    R = Modelica.Constants.R/TmCl_p2.MM,
    G_ref = -801654,
    H_ref = -858975,
    S_ref = -138.072,
    a1 = -2.5363e-06,
    a2 = -3873.88,
    a3 = 0.0003925,
    a4 = -100253,
    c1 = 90.7196,
    c2 = -110052,
    w_ref = 651830,
    z = 2,
    a_0 = 3.72,
    MM = 0.20438) "TmCl+2";

  //1026
  constant DataRecord TmCl2_p1(
    R = Modelica.Constants.R/TmCl2_p1.MM,
    G_ref = -931358,
    H_ref = -1023406,
    S_ref = -78.241,
    a1 = 1.036e-05,
    a2 = -725.59,
    a3 = 0.00026896,
    a4 = -113269,
    c1 = 64.9018,
    c2 = -102717,
    w_ref = 348690,
    z = 1,
    a_0 = 3.72,
    MM = 0.23983) "TmCl2+";

  //1027
  constant DataRecord TmCl3_aq(
    R = Modelica.Constants.R/TmCl3_aq.MM,
    G_ref = -1060226,
    H_ref = -1200808,
    S_ref = -64.015,
    a1 = 2.465e-05,
    a2 = 2762.4,
    a3 = 0.00013213,
    a4 = -127687,
    c1 = -53.1314,
    c2 = -397292,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16893) "TmCl3,aq";

  //1028
  constant DataRecord TmCl4_n1(
    R = Modelica.Constants.R/TmCl4_n1.MM,
    G_ref = -1189093,
    H_ref = -1393690,
    S_ref = -101.253,
    a1 = 4.4796e-05,
    a2 = 7682.24,
    a3 = -6.1404e-05,
    a4 = -148026,
    c1 = -146.686,
    c2 = -993779,
    w_ref = 834750,
    z = -1,
    a_0 = 3.72,
    MM = 0.31073) "TmCl4-";

  //1029
  constant DataRecord TmH2PO4_p2(
    R = Modelica.Constants.R/TmH2PO4_p2.MM,
    G_ref = -1812090,
    H_ref = -2017525,
    S_ref = -163.594,
    a1 = 5.3204e-06,
    a2 = -1953.43,
    a3 = 0.00031665,
    a4 = -108194,
    c1 = 196.1773,
    c2 = 244906,
    w_ref = 688020,
    z = 2,
    a_0 = 3.72,
    MM = 0.26592) "TmH2PO4+2";

  //1030
  constant DataRecord TmNO3_p2(
    R = Modelica.Constants.R/TmNO3_p2.MM,
    G_ref = -781153,
    H_ref = -945584,
    S_ref = -206.271,
    a1 = 4.2376e-06,
    a2 = -2217.98,
    a3 = 0.00032711,
    a4 = -107102,
    c1 = 159.7832,
    c2 = 96215,
    w_ref = 757300,
    z = 2,
    a_0 = 3.72,
    MM = 0.23094) "TmNO3+2";

  //1031
  constant DataRecord TmF_p2(
    R = Modelica.Constants.R/TmF_p2.MM,
    G_ref = -978219,
    H_ref = -1016712,
    S_ref = -84.935,
    a1 = -1.3752e-05,
    a2 = -6611.47,
    a3 = 0.00049993,
    a4 = -88939,
    c1 = 93.3082,
    c2 = -75094,
    w_ref = 570780,
    z = 2,
    a_0 = 3.72,
    MM = 0.18793) "TmF+2";

  //1032
  constant DataRecord TmF2_p1(
    R = Modelica.Constants.R/TmF2_p1.MM,
    G_ref = -1280304,
    H_ref = -1363147,
    S_ref = -66.944,
    a1 = -1.2834e-05,
    a2 = -6387.25,
    a3 = 0.00049111,
    a4 = -89864,
    c1 = 87.2268,
    c2 = -19238,
    w_ref = 330330,
    z = 1,
    a_0 = 3.72,
    MM = 0.20693) "TmF2+";

  //1033
  constant DataRecord TmF3_aq(
    R = Modelica.Constants.R/TmF3_aq.MM,
    G_ref = -1576950,
    H_ref = -1723808,
    S_ref = -115.06,
    a1 = -1.2064e-05,
    a2 = -6199.81,
    a3 = 0.00048384,
    a4 = -90638,
    c1 = -11.252,
    c2 = -251735,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16893) "TmF3,aq";

  //1034
  constant DataRecord TmF4_n1(
    R = Modelica.Constants.R/TmF4_n1.MM,
    G_ref = -1871085,
    H_ref = -2107062,
    S_ref = -246.856,
    a1 = -6.3166e-06,
    a2 = -4794.61,
    a3 = 0.00042823,
    a4 = -96450,
    c1 = -62.9947,
    c2 = -772576,
    w_ref = 1052360,
    z = -1,
    a_0 = 3.72,
    MM = 0.24493) "TmF4-";

  //1035
  constant DataRecord TmOH_p2(
    R = Modelica.Constants.R/TmOH_p2.MM,
    G_ref = -862322,
    H_ref = -913827,
    S_ref = -61.505,
    a1 = 1.0623e-05,
    a2 = -661.28,
    a3 = 0.00026641,
    a4 = -113537,
    c1 = 24.7797,
    c2 = -301683,
    w_ref = 534550,
    z = 2,
    a_0 = 3.72,
    MM = 0.18594) "TmOH+2";

  //1036
  constant DataRecord TmO_p1(
    R = Modelica.Constants.R/TmO_p1.MM,
    G_ref = -815462,
    H_ref = -848013,
    S_ref = 2.092,
    a1 = 1.102e-05,
    a2 = -564.38,
    a3 = 0.00026262,
    a4 = -113935,
    c1 = -103.6937,
    c2 = -650265,
    w_ref = 228610,
    z = 1,
    a_0 = 3.72,
    MM = 0.18493) "TmO+";

  //1037
  constant DataRecord TmO2H_aq(
    R = Modelica.Constants.R/TmO2H_aq.MM,
    G_ref = -1005415,
    H_ref = -1064828,
    S_ref = 145.185,
    a1 = 1.9022e-05,
    a2 = 1388.5,
    a3 = 0.00018609,
    a4 = -122010,
    c1 = -232.3488,
    c2 = -1005666,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.20093) "TmO2H,aq";

  //1038
  constant DataRecord TmO2_n1(
    R = Modelica.Constants.R/TmO2_n1.MM,
    G_ref = -956881,
    H_ref = -1029766,
    S_ref = 99.998,
    a1 = 1.931e-05,
    a2 = 1459.25,
    a3 = 0.00018323,
    a4 = -122303,
    c1 = -130.6044,
    c2 = -840323,
    w_ref = 530070,
    z = -1,
    a_0 = 3.72,
    MM = 0.20093) "TmO2-";

  //1039
  constant DataRecord TmSO4_p1(
    R = Modelica.Constants.R/TmSO4_p1.MM,
    G_ref = -1433857,
    H_ref = -1594606,
    S_ref = -88.701,
    a1 = -5.7501e-06,
    a2 = -1850.5,
    a3 = 0.00031302,
    a4 = -108621,
    c1 = -78.2906,
    c2 = -605094,
    w_ref = 363300,
    z = 1,
    a_0 = 3.72,
    MM = 0.26499) "TmSO4+";

  //1040
  constant DataRecord YbCO3_p1(
    R = Modelica.Constants.R/YbCO3_p1.MM,
    G_ref = -1215452,
    H_ref = -1369883,
    S_ref = -210.037,
    a1 = -4.4396e-06,
    a2 = -4338.6,
    a3 = 0.00041079,
    a4 = -98332,
    c1 = -73.1246,
    c2 = -648558,
    w_ref = 555130,
    z = 1,
    a_0 = 3.72,
    MM = 0.23305) "YbCO3+";

  //1041
  constant DataRecord YbOH_p2(
    R = Modelica.Constants.R/YbOH_p2.MM,
    G_ref = -833871,
    H_ref = -879142,
    S_ref = -54.81,
    a1 = 1.076e-05,
    a2 = -629.27,
    a3 = 0.00026543,
    a4 = -113667,
    c1 = 19.4096,
    c2 = -317022,
    w_ref = 524170,
    z = 2,
    a_0 = 3.72,
    MM = 0.19005) "YbOH+2";

  //1042
  constant DataRecord YbO_p1(
    R = Modelica.Constants.R/YbO_p1.MM,
    G_ref = -787429,
    H_ref = -819144,
    S_ref = 9.205,
    a1 = 1.1153e-05,
    a2 = -531.28,
    a3 = 0.00026117,
    a4 = -114073,
    c1 = -109.4183,
    c2 = -666457,
    w_ref = 217070,
    z = 1,
    a_0 = 3.72,
    MM = 0.18904) "YbO+";

  //1043
  constant DataRecord YbO2H_aq(
    R = Modelica.Constants.R/YbO2H_aq.MM,
    G_ref = -978219,
    H_ref = -1031356,
    S_ref = 152.716,
    a1 = 1.9251e-05,
    a2 = 1445.99,
    a3 = 0.00018351,
    a4 = -122248,
    c1 = -246.587,
    c2 = -1040607,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.20504) "YbO2H,aq";

  //1044
  constant DataRecord YbO2_n1(
    R = Modelica.Constants.R/YbO2_n1.MM,
    G_ref = -928011,
    H_ref = -994311,
    S_ref = 107.947,
    a1 = 1.9557e-05,
    a2 = 1519.8,
    a3 = 0.0001808,
    a4 = -122554,
    c1 = -141.4786,
    c2 = -874414,
    w_ref = 518480,
    z = -1,
    a_0 = 3.72,
    MM = 0.20504) "YbO2-";

  //1045
  constant DataRecord YbHCO3_p2(
    R = Modelica.Constants.R/YbHCO3_p2.MM,
    G_ref = -1237594,
    H_ref = -1355198,
    S_ref = -86.609,
    a1 = 1.439e-07,
    a2 = -3218.46,
    a3 = 0.00036661,
    a4 = -102964,
    c1 = 167.4935,
    c2 = 181544,
    w_ref = 574550,
    z = 2,
    a_0 = 3.72,
    MM = 0.23406) "YbHCO3+2";

  //1046
  constant DataRecord YbCl_p2(
    R = Modelica.Constants.R/YbCl_p2.MM,
    G_ref = -772366,
    H_ref = -823830,
    S_ref = -131.796,
    a1 = -3.5225e-06,
    a2 = -4114.17,
    a3 = 0.00040184,
    a4 = -99261,
    c1 = 78.3337,
    c2 = -150323,
    w_ref = 643160,
    z = 2,
    a_0 = 3.72,
    MM = 0.20849) "YbCl+2";

  //1047
  constant DataRecord YbCl2_p1(
    R = Modelica.Constants.R/YbCl2_p1.MM,
    G_ref = -901234,
    H_ref = -987424,
    S_ref = -70.71,
    a1 = 9.2621e-06,
    a2 = -994.08,
    a3 = 0.00027963,
    a4 = -112160,
    c1 = 41.4358,
    c2 = -181289,
    w_ref = 339360,
    z = 1,
    a_0 = 3.72,
    MM = 0.24394) "YbCl2+";

  //1048
  constant DataRecord YbCl3_aq(
    R = Modelica.Constants.R/YbCl3_aq.MM,
    G_ref = -1029682,
    H_ref = -1163570,
    S_ref = -53.974,
    a1 = 2.3462e-05,
    a2 = 2473.12,
    a3 = 0.00014334,
    a4 = -126495,
    c1 = -53.6832,
    c2 = -530096,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17304) "YbCl3,aq";

  //1049
  constant DataRecord YbCl4_n1(
    R = Modelica.Constants.R/YbCl4_n1.MM,
    G_ref = -1158550,
    H_ref = -1354779,
    S_ref = -87.864,
    a1 = 4.3402e-05,
    a2 = 7342.33,
    a3 = -4.8099e-05,
    a4 = -146624,
    c1 = -206.9829,
    c2 = -1196741,
    w_ref = 814080,
    z = -1,
    a_0 = 3.72,
    MM = 0.31484) "YbCl4-";

  //1050
  constant DataRecord YbH2PO4_p2(
    R = Modelica.Constants.R/YbH2PO4_p2.MM,
    G_ref = -1783639,
    H_ref = -1982798,
    S_ref = -156.105,
    a1 = 4.3321e-06,
    a2 = -2195.76,
    a3 = 0.00032639,
    a4 = -107190,
    c1 = 183.7374,
    c2 = 204639,
    w_ref = 678730,
    z = 2,
    a_0 = 3.72,
    MM = 0.27003) "YbH2PO4+2";

  //1051
  constant DataRecord YbNO3_p2(
    R = Modelica.Constants.R/YbNO3_p2.MM,
    G_ref = -752283,
    H_ref = -910438,
    S_ref = -196.648,
    a1 = 3.2284e-06,
    a2 = -2465.3,
    a3 = 0.00033698,
    a4 = -106077,
    c1 = 146.7638,
    c2 = 55944,
    w_ref = 741740,
    z = 2,
    a_0 = 3.72,
    MM = 0.23505) "YbNO3+2";

  //1052
  constant DataRecord YbF_p2(
    R = Modelica.Constants.R/YbF_p2.MM,
    G_ref = -949350,
    H_ref = -982822,
    S_ref = -81.588,
    a1 = -1.4721e-05,
    a2 = -6848.33,
    a3 = 0.00050926,
    a4 = -87960,
    c1 = 81.3767,
    c2 = -115361,
    w_ref = 567020,
    z = 2,
    a_0 = 3.72,
    MM = 0.19204) "YbF+2";

  //1053
  constant DataRecord YbF2_p1(
    R = Modelica.Constants.R/YbF2_p1.MM,
    G_ref = -1251434,
    H_ref = -1329257,
    S_ref = -63.597,
    a1 = -1.3915e-05,
    a2 = -6649.67,
    a3 = 0.0005011,
    a4 = -88780,
    c1 = 64.2144,
    c2 = -97809,
    w_ref = 325930,
    z = 1,
    a_0 = 3.72,
    MM = 0.21104) "YbF2+";

  //1054
  constant DataRecord YbF3_aq(
    R = Modelica.Constants.R/YbF3_aq.MM,
    G_ref = -1548498,
    H_ref = -1689918,
    S_ref = -110.876,
    a1 = -1.3252e-05,
    a2 = -6489.09,
    a3 = 0.00049504,
    a4 = -89441,
    c1 = -49.4607,
    c2 = -384535,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17304) "YbF3,aq";

  //1055
  constant DataRecord YbF4_n1(
    R = Modelica.Constants.R/YbF4_n1.MM,
    G_ref = -1842634,
    H_ref = -2072335,
    S_ref = -239.743,
    a1 = -7.673e-06,
    a2 = -5125.86,
    a3 = 0.00044125,
    a4 = -95077,
    c1 = -122.2849,
    c2 = -975537,
    w_ref = 1042650,
    z = -1,
    a_0 = 3.72,
    MM = 0.24904) "YbF4-";

  //1056
  constant DataRecord YbSO4_p1(
    R = Modelica.Constants.R/YbSO4_p1.MM,
    G_ref = -1404987,
    H_ref = -1560214,
    S_ref = -84.098,
    a1 = -5.3325e-06,
    a2 = -1952.46,
    a3 = 0.00031703,
    a4 = -108198,
    c1 = -83.8959,
    c2 = -622993,
    w_ref = 358360,
    z = 1,
    a_0 = 3.72,
    MM = 0.2691) "YbSO4+";

  //1057
  constant DataRecord LuCO3_p1(
    R = Modelica.Constants.R/LuCO3_p1.MM,
    G_ref = -1242230,
    H_ref = -1403397,
    S_ref = -241.417,
    a1 = -4.8308e-06,
    a2 = -4431.94,
    a3 = 0.000414,
    a4 = -97947,
    c1 = -64.101,
    c2 = -628956,
    w_ref = 591870,
    z = 1,
    a_0 = 3.72,
    MM = 0.23497) "LuCO3+";

  //1058
  constant DataRecord LuOH_p2(
    R = Modelica.Constants.R/LuOH_p2.MM,
    G_ref = -860649,
    H_ref = -913535,
    S_ref = -89.119,
    a1 = 1.0198e-05,
    a2 = -766.72,
    a3 = 0.00027093,
    a4 = -113098,
    c1 = 47.4495,
    c2 = -236906,
    w_ref = 578350,
    z = 2,
    a_0 = 3.72,
    MM = 0.19197) "LuOH+2";

  //1059
  constant DataRecord LuO_p1(
    R = Modelica.Constants.R/LuO_p1.MM,
    G_ref = -816717,
    H_ref = -851235,
    S_ref = -27.614,
    a1 = 1.0488e-05,
    a2 = -693.79,
    a3 = 0.0002676,
    a4 = -113403,
    c1 = -79.619,
    c2 = -581229,
    w_ref = 274340,
    z = 1,
    a_0 = 3.72,
    MM = 0.19096) "LuO+";

  //1060
  constant DataRecord LuO2H_aq(
    R = Modelica.Constants.R/LuO2H_aq.MM,
    G_ref = -1004997,
    H_ref = -1067338,
    S_ref = 113.386,
    a1 = 1.7934e-05,
    a2 = 1122.61,
    a3 = 0.00019659,
    a4 = -120909,
    c1 = -185.7437,
    c2 = -858222,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.20696) "LuO2H,aq";

  //1061
  constant DataRecord LuO2_n1(
    R = Modelica.Constants.R/LuO2_n1.MM,
    G_ref = -958973,
    H_ref = -1035122,
    S_ref = 66.107,
    a1 = 1.8282e-05,
    a2 = 1209.05,
    a3 = 0.00019288,
    a4 = -121269,
    c1 = -84.9047,
    c2 = -697992,
    w_ref = 581620,
    z = -1,
    a_0 = 3.72,
    MM = 0.20696) "LuO2-";

  //1062
  constant DataRecord LuHCO3_p2(
    R = Modelica.Constants.R/LuHCO3_p2.MM,
    G_ref = -1264823,
    H_ref = -1390762,
    S_ref = -123.428,
    a1 = 1.38e-08,
    a2 = -3249.67,
    a3 = 0.00036771,
    a4 = -102834,
    c1 = 196.9149,
    c2 = 265918,
    w_ref = 630400,
    z = 2,
    a_0 = 3.72,
    MM = 0.23598) "LuHCO3+2";

  //1063
  constant DataRecord LuCl_p2(
    R = Modelica.Constants.R/LuCl_p2.MM,
    G_ref = -797889,
    H_ref = -856046,
    S_ref = -162.758,
    a1 = -3.6903e-06,
    a2 = -4153.75,
    a3 = 0.00040314,
    a4 = -99098,
    c1 = 106.7393,
    c2 = -65944,
    w_ref = 688020,
    z = 2,
    a_0 = 3.72,
    MM = 0.21041) "LuCl+2";

  //1064
  constant DataRecord LuCl2_p1(
    R = Modelica.Constants.R/LuCl2_p1.MM,
    G_ref = -925919,
    H_ref = -1020896,
    S_ref = -109.621,
    a1 = 9.0943e-06,
    a2 = -1036.25,
    a3 = 0.00028152,
    a4 = -111985,
    c1 = 93.9103,
    c2 = -16661,
    w_ref = 394840,
    z = 1,
    a_0 = 3.72,
    MM = 0.24586) "LuCl2+";

  //1065
  constant DataRecord LuCl3_aq(
    R = Modelica.Constants.R/LuCl3_aq.MM,
    G_ref = -1053950,
    H_ref = -1200164,
    S_ref = -105.018,
    a1 = 2.3066e-05,
    a2 = 2375.88,
    a3 = 0.00014729,
    a4 = -126089,
    c1 = -11.2838,
    c2 = -251843,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17496) "LuCl3,aq";

  //1066
  constant DataRecord LuCl4_n1(
    R = Modelica.Constants.R/LuCl4_n1.MM,
    G_ref = -1181980,
    H_ref = -1396619,
    S_ref = -157.737,
    a1 = 4.3319e-05,
    a2 = 7322.38,
    a3 = -4.7375e-05,
    a4 = -146540,
    c1 = -74.8714,
    c2 = -771488,
    w_ref = 920060,
    z = -1,
    a_0 = 3.72,
    MM = 0.31676) "LuCl4-";

  //1067
  constant DataRecord LuH2PO4_p2(
    R = Modelica.Constants.R/LuH2PO4_p2.MM,
    G_ref = -1810835,
    H_ref = -2018362,
    S_ref = -193.719,
    a1 = 4.2091e-06,
    a2 = -2226.18,
    a3 = 0.00032765,
    a4 = -107064,
    c1 = 213.3497,
    c2 = 289014,
    w_ref = 736680,
    z = 2,
    a_0 = 3.72,
    MM = 0.27094) "LuH2PO4+2";

  //1068
  constant DataRecord LuNO3_p2(
    R = Modelica.Constants.R/LuNO3_p2.MM,
    G_ref = -781153,
    H_ref = -951023,
    S_ref = -245.601,
    a1 = 3.1476e-06,
    a2 = -2483.33,
    a3 = 0.00033732,
    a4 = -106002,
    c1 = 177.5313,
    c2 = 140323,
    w_ref = 812240,
    z = 2,
    a_0 = 3.72,
    MM = 0.23697) "LuNO3+2";

  //1069
  constant DataRecord LuF_p2(
    R = Modelica.Constants.R/LuF_p2.MM,
    G_ref = -976127,
    H_ref = -1012110,
    S_ref = -98.742,
    a1 = -1.4949e-05,
    a2 = -6903.35,
    a3 = 0.00051125,
    a4 = -87730,
    c1 = 108.1204,
    c2 = -30987,
    w_ref = 593790,
    z = 2,
    a_0 = 3.72,
    MM = 0.19396) "LuF+2";

  //1070
  constant DataRecord LuF2_p1(
    R = Modelica.Constants.R/LuF2_p1.MM,
    G_ref = -1278630,
    H_ref = -1358963,
    S_ref = -80.333,
    a1 = -1.4177e-05,
    a2 = -6714.19,
    a3 = 0.00050372,
    a4 = -88513,
    c1 = 114.119,
    c2 = 66814,
    w_ref = 353510,
    z = 1,
    a_0 = 3.72,
    MM = 0.21296) "LuF2+";

  //1071
  constant DataRecord LuF3_aq(
    R = Modelica.Constants.R/LuF3_aq.MM,
    G_ref = -1575540,
    H_ref = -1720879,
    S_ref = -133.051,
    a1 = -1.3648e-05,
    a2 = -6586.45,
    a3 = 0.000499,
    a4 = -89040,
    c1 = 30.5942,
    c2 = -106282,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17496) "LuF3,aq";

  //1072
  constant DataRecord LuF4_n1(
    R = Modelica.Constants.R/LuF4_n1.MM,
    G_ref = -1869830,
    H_ref = -2107899,
    S_ref = -276.562,
    a1 = -7.9044e-06,
    a2 = -5182.72,
    a3 = 0.00044356,
    a4 = -94843,
    c1 = 5.7823,
    c2 = -550284,
    w_ref = 1104700,
    z = -1,
    a_0 = 3.72,
    MM = 0.25096) "LuF4-";

  //1073
  constant DataRecord LuSO4_p1(
    R = Modelica.Constants.R/LuSO4_p1.MM,
    G_ref = -1431765,
    H_ref = -1592556,
    S_ref = -111.294,
    a1 = -5.3601e-06,
    a2 = -1945.69,
    a3 = 0.00031676,
    a4 = -108228,
    c1 = -69.2306,
    c2 = -585492,
    w_ref = 400410,
    z = 1,
    a_0 = 3.72,
    MM = 0.27102) "LuSO4+";

  //1074
  constant DataRecord Li_Lac_aq(
    R = Modelica.Constants.R/Li_Lac_aq.MM,
    G_ref = -806407,
    H_ref = -970805,
    S_ref = 128.859,
    a1 = 4.4113e-05,
    a2 = 7516.39,
    a3 = -5.5061e-05,
    a4 = -147344,
    c1 = 330.4327,
    c2 = 935869,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.04801) "Li(Lac),aq";

  //1075
  constant DataRecord Mg_Lac_p1(
    R = Modelica.Constants.R/Mg_Lac_p1.MM,
    G_ref = -974470,
    H_ref = -1148897,
    S_ref = 33.472,
    a1 = 3.1575e-05,
    a2 = 4454.54,
    a3 = 6.5329e-05,
    a4 = -134683,
    c1 = 351.5225,
    c2 = 947249,
    w_ref = 180830,
    z = 1,
    a_0 = 3.72,
    MM = 0.11338) "Mg(Lac)+";

  //1076
  constant DataRecord Mg_Lac2_aq(
    R = Modelica.Constants.R/Mg_Lac2_aq.MM,
    G_ref = -1492960,
    H_ref = -1843889,
    S_ref = 157.879,
    a1 = 7.0939e-05,
    a2 = 14064.43,
    a3 = -0.00031197,
    a4 = -174410,
    c1 = 731.3147,
    c2 = 2329237,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.10645) "Mg(Lac)2,aq";

  //1077
  constant DataRecord Ca_Lac_p1(
    R = Modelica.Constants.R/Ca_Lac_p1.MM,
    G_ref = -1073560,
    H_ref = -1231920,
    S_ref = 96.232,
    a1 = 3.3476e-05,
    a2 = 4917.2,
    a3 = 4.7509e-05,
    a4 = -136595,
    c1 = 330.5841,
    c2 = 905254,
    w_ref = 84730,
    z = 1,
    a_0 = 3.72,
    MM = 0.12915) "Ca(Lac)+";

  //1078
  constant DataRecord Ca_Lac2_aq(
    R = Modelica.Constants.R/Ca_Lac2_aq.MM,
    G_ref = -1592276,
    H_ref = -1921364,
    S_ref = 240.011,
    a1 = 7.342e-05,
    a2 = 14670.32,
    a3 = -0.00033583,
    a4 = -176916,
    c1 = 707.7403,
    c2 = 2247298,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.12222) "Ca(Lac)2,aq";

  //1079
  constant DataRecord Sr_Lac_p1(
    R = Modelica.Constants.R/Sr_Lac_p1.MM,
    G_ref = -1082095,
    H_ref = -1237196,
    S_ref = 121.336,
    a1 = 3.3763e-05,
    a2 = 4988.12,
    a3 = 4.4555e-05,
    a4 = -136892,
    c1 = 313.208,
    c2 = 856929,
    w_ref = 47030,
    z = 1,
    a_0 = 3.72,
    MM = 0.17669) "Sr(Lac)+";

  //1080
  constant DataRecord Sr_Lac2_aq(
    R = Modelica.Constants.R/Sr_Lac2_aq.MM,
    G_ref = -1598698,
    H_ref = -1922217,
    S_ref = 272.86,
    a1 = 7.3882e-05,
    a2 = 14784.25,
    a3 = -0.00034054,
    a4 = -177389,
    c1 = 680.6134,
    c2 = 2153011,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16976) "Sr(Lac)2,aq";

  //1081
  constant DataRecord Ba_Lac_p1(
    R = Modelica.Constants.R/Ba_Lac_p1.MM,
    G_ref = -1077100,
    H_ref = -1219285,
    S_ref = 171.544,
    a1 = 3.6574e-05,
    a2 = 5674.59,
    a3 = 1.7552e-05,
    a4 = -139729,
    c1 = 293.7762,
    c2 = 813780,
    w_ref = -29160,
    z = 1,
    a_0 = 3.72,
    MM = 0.22641) "Ba(Lac)+";

  //1082
  constant DataRecord Ba_Lac2_aq(
    R = Modelica.Constants.R/Ba_Lac2_aq.MM,
    G_ref = -1592104,
    H_ref = -1898088,
    S_ref = 338.565,
    a1 = 7.73e-05,
    a2 = 15619.54,
    a3 = -0.00037349,
    a4 = -180841,
    c1 = 656.3935,
    c2 = 2068825,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.21948) "Ba(Lac)2,aq";

  //1083
  constant DataRecord Mn_Lac_p1(
    R = Modelica.Constants.R/Mn_Lac_p1.MM,
    G_ref = -748857,
    H_ref = -911091,
    S_ref = 73.701,
    a1 = 3.4202e-05,
    a2 = 5094.94,
    a3 = 4.0447e-05,
    a4 = -137331,
    c1 = 352.6204,
    c2 = 971027,
    w_ref = 118490,
    z = 1,
    a_0 = 3.72,
    MM = 0.14401) "Mn(Lac)+";

  //1084
  constant DataRecord Mn_Lac2_aq(
    R = Modelica.Constants.R/Mn_Lac2_aq.MM,
    G_ref = -1267631,
    H_ref = -1602669,
    S_ref = 210.526,
    a1 = 7.4102e-05,
    a2 = 14838.77,
    a3 = -0.00034286,
    a4 = -177611,
    c1 = 744.6633,
    c2 = 2375629,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13708) "Mn(Lac)2,aq";

  //1085
  constant DataRecord Co_Lac_p1(
    R = Modelica.Constants.R/Co_Lac_p1.MM,
    G_ref = -574819,
    H_ref = -752518,
    S_ref = 20.079,
    a1 = 2.9824e-05,
    a2 = 4025.51,
    a3 = 8.2517e-05,
    a4 = -132909,
    c1 = 339.7872,
    c2 = 900075,
    w_ref = 200790,
    z = 1,
    a_0 = 3.72,
    MM = 0.148) "Co(Lac)+";

  //1086
  constant DataRecord Co_Lac2_aq(
    R = Modelica.Constants.R/Co_Lac2_aq.MM,
    G_ref = -1093936,
    H_ref = -1449371,
    S_ref = 140.352,
    a1 = 6.8913e-05,
    a2 = 13571.64,
    a3 = -0.00029301,
    a4 = -172372,
    c1 = 704.8341,
    c2 = 2237193,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14107) "Co(Lac)2,aq";

  //1087
  constant DataRecord Ni_Lac_p1(
    R = Modelica.Constants.R/Ni_Lac_p1.MM,
    G_ref = -567346,
    H_ref = -751534,
    S_ref = -1.598,
    a1 = 2.7004e-05,
    a2 = 3338.83,
    a3 = 0.00010912,
    a4 = -130072,
    c1 = 321.3906,
    c2 = 825286,
    w_ref = 234640,
    z = 1,
    a_0 = 3.72,
    MM = 0.14778) "Ni(Lac)+";

  //1088
  constant DataRecord Ni_Lac2_aq(
    R = Modelica.Constants.R/Ni_Lac2_aq.MM,
    G_ref = -1087660,
    H_ref = -1451413,
    S_ref = 111.985,
    a1 = 6.5644e-05,
    a2 = 12773.08,
    a3 = -0.00026155,
    a4 = -169071,
    c1 = 662.8519,
    c2 = 2091276,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14085) "Ni(Lac)2,aq";

  //1089
  constant DataRecord Cu_Lac_p1(
    R = Modelica.Constants.R/Cu_Lac_p1.MM,
    G_ref = -461692,
    H_ref = -633797,
    S_ref = 41.756,
    a1 = 2.9582e-05,
    a2 = 3968.78,
    a3 = 8.4295e-05,
    a4 = -132675,
    c1 = 348.2674,
    c2 = 940346,
    w_ref = 167070,
    z = 1,
    a_0 = 3.72,
    MM = 0.15261) "Cu(Lac)+";

  //1090
  constant DataRecord Cu_Lac2_aq(
    R = Modelica.Constants.R/Cu_Lac2_aq.MM,
    G_ref = -983437,
    H_ref = -1331282,
    S_ref = 168.72,
    a1 = 6.8771e-05,
    a2 = 13535.83,
    a3 = -0.00029135,
    a4 = -172226,
    c1 = 727.4399,
    c2 = 2315765,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14568) "Cu(Lac)2,aq";

  //1091
  constant DataRecord Zn_Lac_p1(
    R = Modelica.Constants.R/Zn_Lac_p1.MM,
    G_ref = -672499,
    H_ref = -837068,
    S_ref = 75.312,
    a1 = 2.9604e-05,
    a2 = 3973.5,
    a3 = 8.4207e-05,
    a4 = -132696,
    c1 = 340.3266,
    c2 = 928840,
    w_ref = 116820,
    z = 1,
    a_0 = 3.72,
    MM = 0.15444) "Zn(Lac)+";

  //1092
  constant DataRecord Zn_Lac2_aq(
    R = Modelica.Constants.R/Zn_Lac2_aq.MM,
    G_ref = -1194013,
    H_ref = -1526022,
    S_ref = 230.12,
    a1 = 6.8984e-05,
    a2 = 13589.59,
    a3 = -0.00029384,
    a4 = -172448,
    c1 = 720.981,
    c2 = 2293317,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14751) "Zn(Lac)2,aq";

  //1093
  constant DataRecord Cd_Lac_p1(
    R = Modelica.Constants.R/Cd_Lac_p1.MM,
    G_ref = -600023,
    H_ref = -767843,
    S_ref = 74.843,
    a1 = 3.5154e-05,
    a2 = 5328.78,
    a3 = 3.092e-05,
    a4 = -138298,
    c1 = 355.7747,
    c2 = 982533,
    w_ref = 116820,
    z = 1,
    a_0 = 3.72,
    MM = 0.20147) "Cd(Lac)+";

  //1094
  constant DataRecord Cd_Lac2_aq(
    R = Modelica.Constants.R/Cd_Lac2_aq.MM,
    G_ref = -1120053,
    H_ref = -1460572,
    S_ref = 212.016,
    a1 = 7.5168e-05,
    a2 = 15097.29,
    a3 = -0.00035262,
    a4 = -178682,
    c1 = 751.1217,
    c2 = 2398081,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19454) "Cd(Lac)2,aq";

  //1095
  constant DataRecord La_Lac_p2(
    R = Modelica.Constants.R/La_Lac_p2.MM,
    G_ref = -1217908,
    H_ref = -1412941,
    S_ref = -76.726,
    a1 = 2.1978e-05,
    a2 = 2110.03,
    a3 = 0.00015778,
    a4 = -124993,
    c1 = 321.7207,
    c2 = 722376,
    w_ref = 559610,
    z = 2,
    a_0 = 3.72,
    MM = 0.22798) "La(Lac)+2";

  //1096
  constant DataRecord Lu_Lac_p2(
    R = Modelica.Constants.R/Lu_Lac_p2.MM,
    G_ref = -1201569,
    H_ref = -1414531,
    S_ref = -142.591,
    a1 = 1.8237e-05,
    a2 = 1197.38,
    a3 = 0.00019348,
    a4 = -121219,
    c1 = 359.719,
    c2 = 822093,
    w_ref = 660650,
    z = 2,
    a_0 = 3.72,
    MM = 0.26403) "Lu(Lac)+2";

  //1097
  constant DataRecord Na_Lac_aq(
    R = Modelica.Constants.R/Na_Lac_aq.MM,
    G_ref = -774831,
    H_ref = -927463,
    S_ref = 190.51,
    a1 = 4.396e-05,
    a2 = 7478.19,
    a3 = -5.3388e-05,
    a4 = -147185,
    c1 = 302.0739,
    c2 = 837306,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.06406) "Na(Lac),aq";

  //1098
  constant DataRecord Na_Lac2_n1(
    R = Modelica.Constants.R/Na_Lac2_n1.MM,
    G_ref = -1287040,
    H_ref = -1623417,
    S_ref = 290.629,
    a1 = 8.6326e-05,
    a2 = 17821.29,
    a3 = -0.00045962,
    a4 = -189941,
    c1 = 687.5705,
    c2 = 2095845,
    w_ref = 241460,
    z = -1,
    a_0 = 3.72,
    MM = 0.20113) "Na(Lac)2-";

  //1099
  constant DataRecord K_Lac_aq(
    R = Modelica.Constants.R/K_Lac_aq.MM,
    G_ref = -795299,
    H_ref = -935296,
    S_ref = 246.304,
    a1 = 5.0447e-05,
    a2 = 9062.96,
    a3 = -0.0001158,
    a4 = -153737,
    c1 = 263.0125,
    c2 = 701535,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.08017) "K(Lac),aq";

  //1100
  constant DataRecord K_Lac2_n1(
    R = Modelica.Constants.R/K_Lac2_n1.MM,
    G_ref = -1307395,
    H_ref = -1626915,
    S_ref = 360.585,
    a1 = 9.3195e-05,
    a2 = 19499.07,
    a3 = -0.00052562,
    a4 = -196878,
    c1 = 601.5847,
    c2 = 1830948,
    w_ref = 135390,
    z = -1,
    a_0 = 3.72,
    MM = 0.21724) "K(Lac)2-";

  //1101
  constant DataRecord Pb_Lac_p1(
    R = Modelica.Constants.R/Pb_Lac_p1.MM,
    G_ref = -549568,
    H_ref = -684544,
    S_ref = 198.062,
    a1 = 3.4525e-05,
    a2 = 5175.02,
    a3 = 3.6999e-05,
    a4 = -137662,
    c1 = 287.8722,
    c2 = 806110,
    w_ref = -69290,
    z = 1,
    a_0 = 3.72,
    MM = 0.29626) "Pb(Lac)+";

  //1102
  constant DataRecord Pb_Lac2_aq(
    R = Modelica.Constants.R/Pb_Lac2_aq.MM,
    G_ref = -1072338,
    H_ref = -1368670,
    S_ref = 373.267,
    a1 = 7.5168e-05,
    a2 = 15097.29,
    a3 = -0.00035262,
    a4 = -178682,
    c1 = 652.0877,
    c2 = 2053859,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.28933) "Pb(Lac)2,aq";

  //1103
  constant DataRecord Fe_Lac_p1(
    R = Modelica.Constants.R/Fe_Lac_p1.MM,
    G_ref = -616257,
    H_ref = -790830,
    S_ref = 28.066,
    a1 = 3.1184e-05,
    a2 = 4358.39,
    a3 = 6.9316e-05,
    a4 = -134285,
    c1 = 336.9635,
    c2 = 894322,
    w_ref = 188110,
    z = 1,
    a_0 = 3.72,
    MM = 0.14492) "Fe(Lac)+";

  //1104
  constant DataRecord Fe_Lac2_aq(
    R = Modelica.Constants.R/Fe_Lac2_aq.MM,
    G_ref = -1138224,
    H_ref = -1489973,
    S_ref = 150.804,
    a1 = 7.0477e-05,
    a2 = 13953.1,
    a3 = -0.00030791,
    a4 = -173950,
    c1 = 701.6041,
    c2 = 2225972,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13799) "Fe(Lac)2,aq";

  //1105
  constant DataRecord Eu_Lac_p1(
    R = Modelica.Constants.R/Eu_Lac_p1.MM,
    G_ref = -1058129,
    H_ref = -1208352,
    S_ref = 160.41,
    a1 = 4.3266e-05,
    a2 = 7307.69,
    a3 = -4.6497e-05,
    a4 = -146478,
    c1 = 396.293,
    c2 = 1164708,
    w_ref = -12300,
    z = 1,
    a_0 = 3.72,
    MM = 0.24103) "Eu(Lac)+";

  //1106
  constant DataRecord Eu_Lac2_aq(
    R = Modelica.Constants.R/Eu_Lac2_aq.MM,
    G_ref = -1574962,
    H_ref = -1890009,
    S_ref = 323.996,
    a1 = 8.4692e-05,
    a2 = 17424.23,
    a3 = -0.00044437,
    a4 = -188301,
    c1 = 853.3858,
    c2 = 2753524,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.2341) "Eu(Lac)2,aq";

  //1107
  constant DataRecord Li_Glyc_aq(
    R = Modelica.Constants.R/Li_Glyc_aq.MM,
    G_ref = -798947,
    H_ref = -927459,
    S_ref = 112.968,
    a1 = 3.4723e-05,
    a2 = 5222.13,
    a3 = 3.5447e-05,
    a4 = -137859,
    c1 = 238.501,
    c2 = 616341,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.033984) "Li(Glyc),aq";

  //1108
  constant DataRecord Mg_Glyc_p1(
    R = Modelica.Constants.R/Mg_Glyc_p1.MM,
    G_ref = -968550,
    H_ref = -1114827,
    S_ref = -8.368,
    a1 = 2.2398e-05,
    a2 = 2214.05,
    a3 = 0.00015339,
    a4 = -125424,
    c1 = 265.4074,
    c2 = 627717,
    w_ref = 243970,
    z = 1,
    a_0 = 3.72,
    MM = 0.099354) "Mg(Glyc)+";

  //1109
  constant DataRecord Mg_Glyc2_aq(
    R = Modelica.Constants.R/Mg_Glyc2_aq.MM,
    G_ref = -1481065,
    H_ref = -1774183,
    S_ref = 79.278,
    a1 = 5.1087e-05,
    a2 = 9218.98,
    a3 = -0.00012189,
    a4 = -154381,
    c1 = 511.7894,
    c2 = 1566218,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.078397) "Mg(Glyc)2,aq";

  //1110
  constant DataRecord Ca_Glyc_p1(
    R = Modelica.Constants.R/Ca_Glyc_p1.MM,
    G_ref = -1069184,
    H_ref = -1193771,
    S_ref = 73.241,
    a1 = 2.4206e-05,
    a2 = 2654.16,
    a3 = 0.00013637,
    a4 = -127240,
    c1 = 241.9226,
    c2 = 585722,
    w_ref = 120210,
    z = 1,
    a_0 = 3.72,
    MM = 0.11512) "Ca(Glyc)+";

  //1111
  constant DataRecord Ca_Glyc2_aq(
    R = Modelica.Constants.R/Ca_Glyc2_aq.MM,
    G_ref = -1583238,
    H_ref = -1847157,
    S_ref = 186.075,
    a1 = 5.3568e-05,
    a2 = 9824.87,
    a3 = -0.00014577,
    a4 = -156883,
    c1 = 488.2146,
    c2 = 1484278,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.094167) "Ca(Glyc)2,aq";

  //1112
  constant DataRecord Sr_Glyc_p1(
    R = Modelica.Constants.R/Sr_Glyc_p1.MM,
    G_ref = -1078288,
    H_ref = -1196950,
    S_ref = 107.299,
    a1 = 2.4443e-05,
    a2 = 2712.11,
    a3 = 0.00013411,
    a4 = -127482,
    c1 = 223.1829,
    c2 = 537397,
    w_ref = 67700,
    z = 1,
    a_0 = 3.72,
    MM = 0.16266) "Sr(Glyc)+";

  //1113
  constant DataRecord Sr_Glyc2_aq(
    R = Modelica.Constants.R/Sr_Glyc2_aq.MM,
    G_ref = -1590744,
    H_ref = -1845600,
    S_ref = 230.643,
    a1 = 5.403e-05,
    a2 = 9936.21,
    a3 = -0.00014983,
    a4 = -157344,
    c1 = 461.0873,
    c2 = 1389996,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14171) "Sr(Glyc)2,aq";

  //1114
  constant DataRecord Ba_Glyc_p1(
    R = Modelica.Constants.R/Ba_Glyc_p1.MM,
    G_ref = -1073464,
    H_ref = -1183754,
    S_ref = 142.256,
    a1 = 2.7334e-05,
    a2 = 3419.54,
    a3 = 0.00010591,
    a4 = -130407,
    c1 = 205.9101,
    c2 = 494252,
    w_ref = 14980,
    z = 1,
    a_0 = 3.72,
    MM = 0.21238) "Ba(Glyc)+";

  //1115
  constant DataRecord Ba_Glyc2_aq(
    R = Modelica.Constants.R/Ba_Glyc2_aq.MM,
    G_ref = -1584435,
    H_ref = -1827709,
    S_ref = 276.391,
    a1 = 5.7449e-05,
    a2 = 10771.46,
    a3 = -0.00018274,
    a4 = -160799,
    c1 = 436.8669,
    c2 = 1305810,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19143) "Ba(Glyc)2,aq";

  //1116
  constant DataRecord Mn_Glyc_p1(
    R = Modelica.Constants.R/Mn_Glyc_p1.MM,
    G_ref = -744032,
    H_ref = -872757,
    S_ref = 49.852,
    a1 = 2.494e-05,
    a2 = 2835.12,
    a3 = 0.00012891,
    a4 = -127989,
    c1 = 264.1698,
    c2 = 651499,
    w_ref = 156270,
    z = 1,
    a_0 = 3.72,
    MM = 0.12998) "Mn(Glyc)+";

  //1117
  constant DataRecord Mn_Glyc2_aq(
    R = Modelica.Constants.R/Mn_Glyc2_aq.MM,
    G_ref = -1255849,
    H_ref = -1526055,
    S_ref = 155.465,
    a1 = 5.425e-05,
    a2 = 9990.81,
    a3 = -0.00015213,
    a4 = -157569,
    c1 = 525.1372,
    c2 = 1612614,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.10903) "Mn(Glyc)2,aq";

  //1118
  constant DataRecord Co_Glyc_p1(
    R = Modelica.Constants.R/Co_Glyc_p1.MM,
    G_ref = -572668,
    H_ref = -716849,
    S_ref = -3.77,
    a1 = 2.0559e-05,
    a2 = 1763.89,
    a3 = 0.00017131,
    a4 = -123562,
    c1 = 251.258,
    c2 = 580547,
    w_ref = 237690,
    z = 1,
    a_0 = 3.72,
    MM = 0.13397) "Co(Glyc)+";

  //1119
  constant DataRecord Co_Glyc2_aq(
    R = Modelica.Constants.R/Co_Glyc2_aq.MM,
    G_ref = -1088263,
    H_ref = -1378862,
    S_ref = 85.295,
    a1 = 4.9062e-05,
    a2 = 8723.68,
    a3 = -0.00010227,
    a4 = -152331,
    c1 = 485.308,
    c2 = 1474178,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.11302) "Co(Glyc)2,aq";

  //1120
  constant DataRecord Ni_Glyc_p1(
    R = Modelica.Constants.R/Ni_Glyc_p1.MM,
    G_ref = -565480,
    H_ref = -715987,
    S_ref = -25.447,
    a1 = 1.7736e-05,
    a2 = 1075.25,
    a3 = 0.00019824,
    a4 = -120713,
    c1 = 232.7886,
    c2 = 505758,
    w_ref = 270790,
    z = 1,
    a_0 = 3.72,
    MM = 0.13375) "Ni(Glyc)+";

  //1121
  constant DataRecord Ni_Glyc2_aq(
    R = Modelica.Constants.R/Ni_Glyc2_aq.MM,
    G_ref = -1082447,
    H_ref = -1381364,
    S_ref = 56.923,
    a1 = 4.5792e-05,
    a2 = 7925,
    a3 = -7.081e-05,
    a4 = -149030,
    c1 = 443.3253,
    c2 = 1328261,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1128) "Ni(Glyc)2,aq";

  //1122
  constant DataRecord Cu_Glyc_p1(
    R = Modelica.Constants.R/Cu_Glyc_p1.MM,
    G_ref = -457889,
    H_ref = -596475,
    S_ref = 17.908,
    a1 = 2.0315e-05,
    a2 = 1704.98,
    a3 = 0.00017351,
    a4 = -123319,
    c1 = 259.6841,
    c2 = 620814,
    w_ref = 203430,
    z = 1,
    a_0 = 3.72,
    MM = 0.13858) "Cu(Glyc)+";

  //1123
  constant DataRecord Cu_Glyc2_aq(
    R = Modelica.Constants.R/Cu_Glyc2_aq.MM,
    G_ref = -974964,
    H_ref = -1257978,
    S_ref = 113.663,
    a1 = 4.892e-05,
    a2 = 8690.38,
    a3 = -0.00010127,
    a4 = -152197,
    c1 = 507.9138,
    c2 = 1552749,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.11763) "Cu(Glyc)2,aq";

  //1124
  constant DataRecord Zn_Glyc_p1(
    R = Modelica.Constants.R/Zn_Glyc_p1.MM,
    G_ref = -667838,
    H_ref = -814148,
    S_ref = 0.795,
    a1 = 2.0592e-05,
    a2 = 1773.47,
    a3 = 0.00017061,
    a4 = -123600,
    c1 = 258.698,
    c2 = 609308,
    w_ref = 228610,
    z = 1,
    a_0 = 3.72,
    MM = 0.14041) "Zn(Glyc)+";

  //1125
  constant DataRecord Zn_Glyc2_aq(
    R = Modelica.Constants.R/Zn_Glyc2_aq.MM,
    G_ref = -1185373,
    H_ref = -1477534,
    S_ref = 91.266,
    a1 = 4.9133e-05,
    a2 = 8741.63,
    a3 = -0.00010311,
    a4 = -152406,
    c1 = 501.4545,
    c2 = 1530302,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.11946) "Zn(Glyc)2,aq";

  //1126
  constant DataRecord Cd_Glyc_p1(
    R = Modelica.Constants.R/Cd_Glyc_p1.MM,
    G_ref = -595304,
    H_ref = -729610,
    S_ref = 50.995,
    a1 = 2.5891e-05,
    a2 = 3066.58,
    a3 = 0.0001199,
    a4 = -128947,
    c1 = 267.289,
    c2 = 663001,
    w_ref = 154220,
    z = 1,
    a_0 = 3.72,
    MM = 0.18744) "Cd(Glyc)+";

  //1127
  constant DataRecord Cd_Glyc2_aq(
    R = Modelica.Constants.R/Cd_Glyc2_aq.MM,
    G_ref = -1110383,
    H_ref = -1386071,
    S_ref = 156.959,
    a1 = 5.5317e-05,
    a2 = 10251.8,
    a3 = -0.00016255,
    a4 = -158649,
    c1 = 531.5964,
    c2 = 1635061,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16649) "Cd(Glyc)2,aq";

  //1128
  constant DataRecord Eu_Glyc_p1(
    R = Modelica.Constants.R/Eu_Glyc_p1.MM,
    G_ref = -1054552,
    H_ref = -1171261,
    S_ref = 136.562,
    a1 = 3.3997e-05,
    a2 = 5044.31,
    a3 = 4.2535e-05,
    a4 = -137122,
    c1 = 307.6688,
    c2 = 845176,
    w_ref = 23600,
    z = 1,
    a_0 = 3.72,
    MM = 0.227) "Eu(Glyc)+";

  //1129
  constant DataRecord Eu_Glyc2_aq(
    R = Modelica.Constants.R/Eu_Glyc2_aq.MM,
    G_ref = -1565008,
    H_ref = -1815224,
    S_ref = 268.939,
    a1 = 6.4841e-05,
    a2 = 12576.14,
    a3 = -0.00025363,
    a4 = -168260,
    c1 = 633.8597,
    c2 = 1990509,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.20605) "Eu(Glyc)2,aq";

  //1130
  constant DataRecord Na_Glyc_aq(
    R = Modelica.Constants.R/Na_Glyc_aq.MM,
    G_ref = -769141,
    H_ref = -888259,
    S_ref = 166.661,
    a1 = 3.457e-05,
    a2 = 5186.57,
    a3 = 3.6459e-05,
    a4 = -137712,
    c1 = 210.1427,
    c2 = 517774,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.050034) "Na(Glyc),aq";

  //1131
  constant DataRecord Na_Glyc2_n1(
    R = Modelica.Constants.R/Na_Glyc2_n1.MM,
    G_ref = -1275601,
    H_ref = -1546758,
    S_ref = 236.877,
    a1 = 6.6749e-05,
    a2 = 13042.28,
    a3 = -0.00027199,
    a4 = -170184,
    c1 = 475.538,
    c2 = 1332830,
    w_ref = 322840,
    z = -1,
    a_0 = 3.72,
    MM = 0.17308) "Na(Glyc)2-";

  //1132
  constant DataRecord K_Glyc_aq(
    R = Modelica.Constants.R/K_Glyc_aq.MM,
    G_ref = -789609,
    H_ref = -896091,
    S_ref = 222.455,
    a1 = 4.1058e-05,
    a2 = 6768.62,
    a3 = -2.5288e-05,
    a4 = -144252,
    c1 = 171.0804,
    c2 = 382008,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.066144) "K(Glyc),aq";

  //1133
  constant DataRecord K_Glyc2_n1(
    R = Modelica.Constants.R/K_Glyc2_n1.MM,
    G_ref = -1295956,
    H_ref = -1550251,
    S_ref = 306.834,
    a1 = 7.362e-05,
    a2 = 14719.73,
    a3 = -0.0003379,
    a4 = -177121,
    c1 = 389.5689,
    c2 = 1067933,
    w_ref = 216940,
    z = -1,
    a_0 = 3.72,
    MM = 0.18919) "K(Glyc)2-";

  //1134
  constant DataRecord Pb_Glyc_p1(
    R = Modelica.Constants.R/Pb_Glyc_p1.MM,
    G_ref = -543995,
    H_ref = -645453,
    S_ref = 174.213,
    a1 = 3.4647e-05,
    a2 = 5205.57,
    a3 = 3.5685e-05,
    a4 = -137787,
    c1 = 199.2864,
    c2 = 486582,
    w_ref = -32970,
    z = 1,
    a_0 = 3.72,
    MM = 0.28223) "Pb(Glyc)+";

  //1135
  constant DataRecord Pb_Glyc2_aq(
    R = Modelica.Constants.R/Pb_Glyc2_aq.MM,
    G_ref = -1061129,
    H_ref = -1292630,
    S_ref = 318.423,
    a1 = 7.5168e-05,
    a2 = 15097.29,
    a3 = -0.00035262,
    a4 = -178682,
    c1 = 432.5612,
    c2 = 1290843,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.26128) "Pb(Glyc)2,aq";

  //1136
  constant DataRecord Fe_Glyc_p1(
    R = Modelica.Constants.R/Fe_Glyc_p1.MM,
    G_ref = -610906,
    H_ref = -751915,
    S_ref = 4.217,
    a1 = 3.1311e-05,
    a2 = 4390.35,
    a3 = 6.7852e-05,
    a4 = -134419,
    c1 = 248.4945,
    c2 = 574794,
    w_ref = 225680,
    z = 1,
    a_0 = 3.72,
    MM = 0.13089) "Fe(Glyc)+";

  //1137
  constant DataRecord Fe_Glyc2_aq(
    R = Modelica.Constants.R/Fe_Glyc2_aq.MM,
    G_ref = -1127416,
    H_ref = -1414234,
    S_ref = 95.742,
    a1 = 7.0477e-05,
    a2 = 13953.1,
    a3 = -0.00030791,
    a4 = -173950,
    c1 = 482.0788,
    c2 = 1462952,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.10994) "Fe(Glyc)2,aq";

  //1138
  constant DataRecord Cd_Alan_p1(
    R = Modelica.Constants.R/Cd_Alan_p1.MM,
    G_ref = -418743,
    H_ref = -590011,
    S_ref = 122.332,
    a1 = 3.8138e-05,
    a2 = 6056.01,
    a3 = 2.636e-06,
    a4 = -141306,
    c1 = 308.7432,
    c2 = 841750,
    w_ref = 45940,
    z = 1,
    a_0 = 3.72,
    MM = 0.20048) "Cd(Alan)+";

  //1139
  constant DataRecord Cd_Alan2_aq(
    R = Modelica.Constants.R/Cd_Alan2_aq.MM,
    G_ref = -752869,
    H_ref = -1102149,
    S_ref = 300.733,
    a1 = 8.1983e-05,
    a2 = 16763.32,
    a3 = -0.00041854,
    a4 = -185569,
    c1 = 654.3998,
    c2 = 2061900,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22457) "Cd(Alan)2,aq";

  //1140
  constant DataRecord Ca_Alan_p1(
    R = Modelica.Constants.R/Ca_Alan_p1.MM,
    G_ref = -872247,
    H_ref = -1033795,
    S_ref = 144.578,
    a1 = 3.6453e-05,
    a2 = 5643.71,
    a3 = 1.8995e-05,
    a4 = -139599,
    c1 = 283.3518,
    c2 = 764471,
    w_ref = 11670,
    z = 1,
    a_0 = 3.72,
    MM = 0.12816) "Ca(Alan)+";

  //1141
  constant DataRecord Ca_Alan2_aq(
    R = Modelica.Constants.R/Ca_Alan2_aq.MM,
    G_ref = -1188448,
    H_ref = -1525963,
    S_ref = 329.846,
    a1 = 8.0234e-05,
    a2 = 16336.34,
    a3 = -0.00040175,
    a4 = -183803,
    c1 = 611.0184,
    c2 = 1911117,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15225) "Ca(Alan)2,aq";

  //1142
  constant DataRecord Pb_Alan_p1(
    R = Modelica.Constants.R/Pb_Alan_p1.MM,
    G_ref = -364807,
    H_ref = -503231,
    S_ref = 245.551,
    a1 = 3.7505e-05,
    a2 = 5903.12,
    a3 = 8.335e-06,
    a4 = -140674,
    c1 = 240.7528,
    c2 = 665327,
    w_ref = -141080,
    z = 1,
    a_0 = 3.72,
    MM = 0.29527) "Pb(Alan)+";

  //1143
  constant DataRecord Pb_Alan2_aq(
    R = Modelica.Constants.R/Pb_Alan2_aq.MM,
    G_ref = -695678,
    H_ref = -1000775,
    S_ref = 461.981,
    a1 = 8.1983e-05,
    a2 = 16763.32,
    a3 = -0.00041854,
    a4 = -185569,
    c1 = 555.3658,
    c2 = 1717678,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.31936) "Pb(Alan)2,aq";

  //1144
  constant DataRecord Mg_Alan_p1(
    R = Modelica.Constants.R/Mg_Alan_p1.MM,
    G_ref = -777550,
    H_ref = -969621,
    S_ref = 33.338,
    a1 = 3.4798e-05,
    a2 = 5241.8,
    a3 = 3.4351e-05,
    a4 = -137938,
    c1 = 311.0181,
    c2 = 806466,
    w_ref = 180830,
    z = 1,
    a_0 = 3.72,
    MM = 0.11239) "Mg(Alan)+";

  //1145
  constant DataRecord Mg_Alan2_aq(
    R = Modelica.Constants.R/Mg_Alan2_aq.MM,
    G_ref = -1097175,
    H_ref = -1475450,
    S_ref = 184.272,
    a1 = 7.7754e-05,
    a2 = 15730.46,
    a3 = -0.00037789,
    a4 = -181297,
    c1 = 634.5927,
    c2 = 1993057,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13648) "Mg(Alan)2,aq";

  //1146
  constant DataRecord Sr_Alan_p1(
    R = Modelica.Constants.R/Sr_Alan_p1.MM,
    G_ref = -880439,
    H_ref = -1036059,
    S_ref = 178.636,
    a1 = 3.6694e-05,
    a2 = 5703.29,
    a3 = 1.6514e-05,
    a4 = -139846,
    c1 = 264.7179,
    c2 = 716146,
    w_ref = -39660,
    z = 1,
    a_0 = 3.72,
    MM = 0.1757) "Sr(Alan)+";

  //1147
  constant DataRecord Sr_Alan2_aq(
    R = Modelica.Constants.R/Sr_Alan2_aq.MM,
    G_ref = -1194243,
    H_ref = -1522696,
    S_ref = 374.414,
    a1 = 8.0696e-05,
    a2 = 16447.68,
    a3 = -0.00040581,
    a4 = -184263,
    c1 = 583.8914,
    c2 = 1816831,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19979) "Sr(Alan)2,aq";

  //1148
  constant DataRecord Mn_Alan_p1(
    R = Modelica.Constants.R/Mn_Alan_p1.MM,
    G_ref = -558899,
    H_ref = -724585,
    S_ref = 121.19,
    a1 = 3.7184e-05,
    a2 = 5822.75,
    a3 = 1.1891e-05,
    a4 = -140340,
    c1 = 305.5312,
    c2 = 830244,
    w_ref = 47030,
    z = 1,
    a_0 = 3.72,
    MM = 0.14302) "Mn(Alan)+";

  //1149
  constant DataRecord Mn_Alan2_aq(
    R = Modelica.Constants.R/Mn_Alan2_aq.MM,
    G_ref = -887318,
    H_ref = -1231121,
    S_ref = 299.24,
    a1 = 8.0916e-05,
    a2 = 16502.28,
    a3 = -0.0004081,
    a4 = -184489,
    c1 = 647.9414,
    c2 = 2039449,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16711) "Mn(Alan)2,aq";

  //1150
  constant DataRecord Co_Alan_p1(
    R = Modelica.Constants.R/Co_Alan_p1.MM,
    G_ref = -393710,
    H_ref = -570049,
    S_ref = 83.68,
    a1 = 3.272e-05,
    a2 = 4734.11,
    a3 = 5.435e-05,
    a4 = -135842,
    c1 = 290.3508,
    c2 = 759291,
    w_ref = 103850,
    z = 1,
    a_0 = 3.72,
    MM = 0.14701) "Co(Alan)+";

  //1151
  constant DataRecord Co_Alan2_aq(
    R = Modelica.Constants.R/Co_Alan2_aq.MM,
    G_ref = -727150,
    H_ref = -1084794,
    S_ref = 251.04,
    a1 = 7.5728e-05,
    a2 = 15235.07,
    a3 = -0.00035826,
    a4 = -179251,
    c1 = 608.1122,
    c2 = 1901013,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1711) "Co(Alan)2,aq";

  //1152
  constant DataRecord Ni_Alan_p1(
    R = Modelica.Constants.R/Ni_Alan_p1.MM,
    G_ref = -391317,
    H_ref = -573756,
    S_ref = 62.76,
    a1 = 2.9895e-05,
    a2 = 4043.42,
    a3 = 8.1739e-05,
    a4 = -132984,
    c1 = 271.8357,
    c2 = 684502,
    w_ref = 136400,
    z = 1,
    a_0 = 3.72,
    MM = 0.14679) "Ni(Alan)+";

  //1153
  constant DataRecord Ni_Alan2_aq(
    R = Modelica.Constants.R/Ni_Alan2_aq.MM,
    G_ref = -730296,
    H_ref = -1100275,
    S_ref = 209.2,
    a1 = 7.2458e-05,
    a2 = 14436.47,
    a3 = -0.00032679,
    a4 = -175950,
    c1 = 566.1299,
    c2 = 1755096,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17088) "Ni(Alan)2,aq";

  //1154
  constant DataRecord Cu_Alan_p1(
    R = Modelica.Constants.R/Cu_Alan_p1.MM,
    G_ref = -295369,
    H_ref = -460114,
    S_ref = 104.6,
    a1 = 3.2379e-05,
    a2 = 4651.14,
    a3 = 5.758e-05,
    a4 = -135499,
    c1 = 296.1264,
    c2 = 799562,
    w_ref = 40750,
    z = 1,
    a_0 = 3.72,
    MM = 0.15162) "Cu(Alan)+";

  //1155
  constant DataRecord Cu_Alan2_aq(
    R = Modelica.Constants.R/Cu_Alan2_aq.MM,
    G_ref = -647072,
    H_ref = -993114,
    S_ref = 271.96,
    a1 = 7.5586e-05,
    a2 = 15199.26,
    a3 = -0.0003566,
    a4 = -179104,
    c1 = 630.7179,
    c2 = 1979584,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17571) "Cu(Alan)2,aq";

  //1156
  constant DataRecord Zn_Alan_p1(
    R = Modelica.Constants.R/Zn_Alan_p1.MM,
    G_ref = -487909,
    H_ref = -671381,
    S_ref = 71.128,
    a1 = 3.285e-05,
    a2 = 4765.2,
    a3 = 5.3329e-05,
    a4 = -135967,
    c1 = 300.4556,
    c2 = 788056,
    w_ref = 123680,
    z = 1,
    a_0 = 3.72,
    MM = 0.15345) "Zn(Alan)+";

  //1157
  constant DataRecord Zn_Alan2_aq(
    R = Modelica.Constants.R/Zn_Alan2_aq.MM,
    G_ref = -824717,
    H_ref = -1185700,
    S_ref = 251.04,
    a1 = 7.5799e-05,
    a2 = 15253.11,
    a3 = -0.00035908,
    a4 = -179326,
    c1 = 624.2591,
    c2 = 1957137,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17754) "Zn(Alan)2,aq";

  //1158
  constant DataRecord Ba_Alan_p1(
    R = Modelica.Constants.R/Ba_Alan_p1.MM,
    G_ref = -878694,
    H_ref = -1019653,
    S_ref = 234.71,
    a1 = 3.9475e-05,
    a2 = 6381.77,
    a3 = -1.0004e-05,
    a4 = -142653,
    c1 = 244.4816,
    c2 = 672996,
    w_ref = -124560,
    z = 1,
    a_0 = 3.72,
    MM = 0.22542) "Ba(Alan)+";

  //1159
  constant DataRecord Ba_Alan2_aq(
    R = Modelica.Constants.R/Ba_Alan2_aq.MM,
    G_ref = -1193641,
    H_ref = -1502269,
    S_ref = 447.797,
    a1 = 8.4115e-05,
    a2 = 17282.97,
    a3 = -0.00043876,
    a4 = -187715,
    c1 = 559.6715,
    c2 = 1732645,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.24951) "Ba(Alan)2,aq";

  //1160
  constant DataRecord Fe_Alan_p1(
    R = Modelica.Constants.R/Fe_Alan_p1.MM,
    G_ref = -432006,
    H_ref = -609885,
    S_ref = 75.555,
    a1 = 3.4167e-05,
    a2 = 5085.99,
    a3 = 4.0836e-05,
    a4 = -137294,
    c1 = 289.891,
    c2 = 753538,
    w_ref = 116820,
    z = 1,
    a_0 = 3.72,
    MM = 0.14393) "Fe(Alan)+";

  //1161
  constant DataRecord Fe_Alan2_aq(
    R = Modelica.Constants.R/Fe_Alan2_aq.MM,
    G_ref = -765618,
    H_ref = -1123550,
    S_ref = 239.517,
    a1 = 7.7292e-05,
    a2 = 15616.53,
    a3 = -0.00037317,
    a4 = -180828,
    c1 = 604.8821,
    c2 = 1889791,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16802) "Fe(Alan)2,aq";

  //1162
  constant DataRecord Eu_Alan_p1(
    R = Modelica.Constants.R/Eu_Alan_p1.MM,
    G_ref = -859097,
    H_ref = -1012766,
    S_ref = 207.899,
    a1 = 4.6247e-05,
    a2 = 8035.71,
    a3 = -7.509e-05,
    a4 = -149490,
    c1 = 349.1916,
    c2 = 1023925,
    w_ref = -83930,
    z = 1,
    a_0 = 3.72,
    MM = 0.24004) "Eu(Alan)+";

  //1163
  constant DataRecord Eu_Alan2_aq(
    R = Modelica.Constants.R/Eu_Alan2_aq.MM,
    G_ref = -1176210,
    H_ref = -1500023,
    S_ref = 412.71,
    a1 = 9.1507e-05,
    a2 = 19087.62,
    a3 = -0.00050961,
    a4 = -195179,
    c1 = 756.6638,
    c2 = 2417344,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.26413) "Eu(Alan)2,aq";

  //1164
  constant DataRecord Cd_Gly_p1(
    R = Modelica.Constants.R/Cd_Gly_p1.MM,
    G_ref = -419387,
    H_ref = -552656,
    S_ref = 112.968,
    a1 = 2.7814e-05,
    a2 = 3536.53,
    a3 = 0.00010138,
    a4 = -130888,
    c1 = 198.9266,
    c2 = 455784,
    w_ref = 59290,
    z = 1,
    a_0 = 3.72,
    MM = 0.18646) "Cd(Gly)+";

  //1165
  constant DataRecord Cd_Gly2_aq(
    R = Modelica.Constants.R/Cd_Gly2_aq.MM,
    G_ref = -755530,
    H_ref = -1031804,
    S_ref = 271.96,
    a1 = 6.0061e-05,
    a2 = 11410.69,
    a3 = -0.00020814,
    a4 = -163440,
    c1 = 389.23,
    c2 = 1140236,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19652) "Cd(Gly)2,aq";

  //1166
  constant DataRecord Ca_Gly_p1(
    R = Modelica.Constants.R/Ca_Gly_p1.MM,
    G_ref = -875460,
    H_ref = -998424,
    S_ref = 137.172,
    a1 = 2.6121e-05,
    a2 = 3121.05,
    a3 = 0.00011816,
    a4 = -129173,
    c1 = 173.3205,
    c2 = 378501,
    w_ref = 22720,
    z = 1,
    a_0 = 3.72,
    MM = 0.11414) "Ca(Gly)+";

  //1167
  constant DataRecord Ca_Gly2_aq(
    R = Modelica.Constants.R/Ca_Gly2_aq.MM,
    G_ref = -1194762,
    H_ref = -1455789,
    S_ref = 312.75,
    a1 = 5.8313e-05,
    a2 = 10983.71,
    a3 = -0.00019135,
    a4 = -161674,
    c1 = 345.8482,
    c2 = 989453,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1242) "Ca(Gly)2,aq";

  //1168
  constant DataRecord Sr_Gly_p1(
    R = Modelica.Constants.R/Sr_Gly_p1.MM,
    G_ref = -884221,
    H_ref = -1001260,
    S_ref = 171.23,
    a1 = 2.6362e-05,
    a2 = 3180.63,
    a3 = 0.0001157,
    a4 = -129419,
    c1 = 154.6904,
    c2 = 330176,
    w_ref = -28620,
    z = 1,
    a_0 = 3.72,
    MM = 0.16168) "Sr(Gly)+";

  //1169
  constant DataRecord Sr_Gly2_aq(
    R = Modelica.Constants.R/Sr_Gly2_aq.MM,
    G_ref = -1201641,
    H_ref = -1453605,
    S_ref = 357.318,
    a1 = 5.8775e-05,
    a2 = 11094.96,
    a3 = -0.00019541,
    a4 = -162134,
    c1 = 318.7216,
    c2 = 895167,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17174) "Sr(Gly)2,aq";

  //1170
  constant DataRecord Mn_Gly_p1(
    R = Modelica.Constants.R/Mn_Gly_p1.MM,
    G_ref = -563882,
    H_ref = -693720,
    S_ref = 104.6,
    a1 = 2.6903e-05,
    a2 = 3313.14,
    a3 = 0.00011034,
    a4 = -129968,
    c1 = 196.8543,
    c2 = 444278,
    w_ref = 72760,
    z = 1,
    a_0 = 3.72,
    MM = 0.129) "Mn(Gly)+";

  //1171
  constant DataRecord Mn_Gly2_aq(
    R = Modelica.Constants.R/Mn_Gly2_aq.MM,
    G_ref = -895799,
    H_ref = -1166696,
    S_ref = 270.123,
    a1 = 5.8995e-05,
    a2 = 11149.57,
    a3 = -0.00019771,
    a4 = -162360,
    c1 = 382.7716,
    c2 = 1117785,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13906) "Mn(Gly)2,aq";

  //1172
  constant DataRecord Fe_Gly_p1(
    R = Modelica.Constants.R/Fe_Gly_p1.MM,
    G_ref = -431111,
    H_ref = -565890,
    S_ref = 83.68,
    a1 = 2.3754e-05,
    a2 = 2545.42,
    a3 = 0.0001403,
    a4 = -126792,
    c1 = 177.6497,
    c2 = 367573,
    w_ref = 103850,
    z = 1,
    a_0 = 3.72,
    MM = 0.12991) "Fe(Gly)+";

  //1173
  constant DataRecord Fe_Gly2_aq(
    R = Modelica.Constants.R/Fe_Gly2_aq.MM,
    G_ref = -764454,
    H_ref = -1039837,
    S_ref = 242.747,
    a1 = 5.5371e-05,
    a2 = 10263.81,
    a3 = -0.00016277,
    a4 = -158699,
    c1 = 339.7123,
    c2 = 968127,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13997) "Fe(Gly)2,aq";

  //1174
  constant DataRecord Co_Gly_p1(
    R = Modelica.Constants.R/Co_Gly_p1.MM,
    G_ref = -398296,
    H_ref = -540079,
    S_ref = 62.76,
    a1 = 2.2461e-05,
    a2 = 2228.94,
    a3 = 0.00015288,
    a4 = -125482,
    c1 = 182.3061,
    c2 = 373326,
    w_ref = 136400,
    z = 1,
    a_0 = 3.72,
    MM = 0.13299) "Co(Gly)+";

  //1175
  constant DataRecord Co_Gly2_aq(
    R = Modelica.Constants.R/Co_Gly2_aq.MM,
    G_ref = -736204,
    H_ref = -1018499,
    S_ref = 230.12,
    a1 = 5.3807e-05,
    a2 = 9882.48,
    a3 = -0.00014787,
    a4 = -157122,
    c1 = 342.9424,
    c2 = 979349,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14305) "Co(Gly)2,aq";

  //1176
  constant DataRecord Pb_Gly_p1(
    R = Modelica.Constants.R/Pb_Gly_p1.MM,
    G_ref = -370075,
    H_ref = -469913,
    S_ref = 238.145,
    a1 = 2.7175e-05,
    a2 = 3380.25,
    a3 = 0.00010758,
    a4 = -130244,
    c1 = 130.7429,
    c2 = 279361,
    w_ref = -129830,
    z = 1,
    a_0 = 3.72,
    MM = 0.28125) "Pb(Gly)+";

  //1177
  constant DataRecord Pb_Gly2_aq(
    R = Modelica.Constants.R/Pb_Gly2_aq.MM,
    G_ref = -704389,
    H_ref = -932999,
    S_ref = 444.885,
    a1 = 6.0061e-05,
    a2 = 11410.69,
    a3 = -0.00020814,
    a4 = -163440,
    c1 = 290.1947,
    c2 = 796019,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.29131) "Pb(Gly)2,aq";

  //1178
  constant DataRecord Mg_Gly_p1(
    R = Modelica.Constants.R/Mg_Gly_p1.MM,
    G_ref = -788642,
    H_ref = -942128,
    S_ref = 25.932,
    a1 = 2.4462e-05,
    a2 = 2717.63,
    a3 = 0.00013365,
    a4 = -127503,
    c1 = 200.868,
    c2 = 420500,
    w_ref = 190580,
    z = 1,
    a_0 = 3.72,
    MM = 0.098368) "Mg(Gly)+";

  //1179
  constant DataRecord Mg_Gly2_aq(
    R = Modelica.Constants.R/Mg_Gly2_aq.MM,
    G_ref = -1120785,
    H_ref = -1422573,
    S_ref = 167.176,
    a1 = 5.5833e-05,
    a2 = 10377.74,
    a3 = -0.00016748,
    a4 = -159172,
    c1 = 369.4229,
    c2 = 1071393,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.10843) "Mg(Gly)2,aq";

  //1180
  constant DataRecord Ni_Gly_p1(
    R = Modelica.Constants.R/Ni_Gly_p1.MM,
    G_ref = -395560,
    H_ref = -540945,
    S_ref = 50.208,
    a1 = 1.9587e-05,
    a2 = 1527.79,
    a3 = 0.00018032,
    a4 = -122587,
    c1 = 162.4287,
    c2 = 298537,
    w_ref = 154220,
    z = 1,
    a_0 = 3.72,
    MM = 0.13277) "Ni(Gly)+";

  //1181
  constant DataRecord Ni_Gly2_aq(
    R = Modelica.Constants.R/Ni_Gly2_aq.MM,
    G_ref = -738606,
    H_ref = -1029494,
    S_ref = 200.832,
    a1 = 5.0537e-05,
    a2 = 9083.76,
    a3 = -0.0001164,
    a4 = -153821,
    c1 = 300.9597,
    c2 = 833432,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14283) "Ni(Gly)2,aq";

  //1182
  constant DataRecord Cu_Gly_p1(
    R = Modelica.Constants.R/Cu_Gly_p1.MM,
    G_ref = -298298,
    H_ref = -428475,
    S_ref = 104.6,
    a1 = 2.2118e-05,
    a2 = 2143.8,
    a3 = 0.00015654,
    a4 = -125131,
    c1 = 188.0269,
    c2 = 413597,
    w_ref = 72760,
    z = 1,
    a_0 = 3.72,
    MM = 0.1376) "Cu(Gly)+";

  //1183
  constant DataRecord Cu_Gly2_aq(
    R = Modelica.Constants.R/Cu_Gly2_aq.MM,
    G_ref = -654700,
    H_ref = -927886,
    S_ref = 263.592,
    a1 = 5.3665e-05,
    a2 = 9846.54,
    a3 = -0.00014621,
    a4 = -156975,
    c1 = 365.5473,
    c2 = 1057924,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14766) "Cu(Gly)2,aq";

  //1184
  constant DataRecord Zn_Gly_p1(
    R = Modelica.Constants.R/Zn_Gly_p1.MM,
    G_ref = -492951,
    H_ref = -634332,
    S_ref = 75.312,
    a1 = 2.2458e-05,
    a2 = 2226.98,
    a3 = 0.00015325,
    a4 = -125474,
    c1 = 188.7758,
    c2 = 402091,
    w_ref = 116820,
    z = 1,
    a_0 = 3.72,
    MM = 0.13943) "Zn(Gly)+";

  //1185
  constant DataRecord Zn_Gly2_aq(
    R = Modelica.Constants.R/Zn_Gly2_aq.MM,
    G_ref = -833198,
    H_ref = -1118835,
    S_ref = 230.12,
    a1 = 5.3878e-05,
    a2 = 9900.39,
    a3 = -0.00014869,
    a4 = -157197,
    c1 = 359.0889,
    c2 = 1035473,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.12547) "Zn(Gly)2,aq";

  //1186
  constant DataRecord Eu_Gly_p1(
    R = Modelica.Constants.R/Eu_Gly_p1.MM,
    G_ref = -864536,
    H_ref = -979625,
    S_ref = 200.493,
    a1 = 3.5916e-05,
    a2 = 5512.84,
    a3 = 2.4129e-05,
    a4 = -139059,
    c1 = 239.1763,
    c2 = 637960,
    w_ref = -72720,
    z = 1,
    a_0 = 3.72,
    MM = 0.22602) "Eu(Gly)+";

  //1187
  constant DataRecord Eu_Gly2_aq(
    R = Modelica.Constants.R/Eu_Gly2_aq.MM,
    G_ref = -1187490,
    H_ref = -1434815,
    S_ref = 395.614,
    a1 = 6.9586e-05,
    a2 = 13734.9,
    a3 = -0.00029922,
    a4 = -173050,
    c1 = 491.4941,
    c2 = 1495680,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.23608) "Eu(Gly)2,aq";

  //1188
  constant DataRecord Ba_Gly_p1(
    R = Modelica.Constants.R/Ba_Gly_p1.MM,
    G_ref = -884251,
    H_ref = -986621,
    S_ref = 227.304,
    a1 = 2.9144e-05,
    a2 = 3861.62,
    a3 = 8.855e-05,
    a4 = -132231,
    c1 = 134.4692,
    c2 = 287031,
    w_ref = -113340,
    z = 1,
    a_0 = 3.72,
    MM = 0.2114) "Ba(Gly)+";

  //1189
  constant DataRecord Ba_Gly2_aq(
    R = Modelica.Constants.R/Ba_Gly2_aq.MM,
    G_ref = -1204235,
    H_ref = -1436376,
    S_ref = 430.701,
    a1 = 6.2194e-05,
    a2 = 11930.26,
    a3 = -0.00022834,
    a4 = -165590,
    c1 = 294.5009,
    c2 = 810985,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22146) "Ba(Gly)2,aq";

  //1190
  constant DataRecord Mg_For_p1(
    R = Modelica.Constants.R/Mg_For_p1.MM,
    G_ref = -813027,
    H_ref = -902397,
    S_ref = -56.827,
    a1 = 1.4808e-05,
    a2 = 359.53,
    a3 = 0.00022652,
    a4 = -117754,
    c1 = 142.4585,
    c2 = 176895,
    w_ref = 317310,
    z = 1,
    a_0 = 3.72,
    MM = 0.069328) "Mg(For)+";

  //1191
  constant DataRecord Mg_For2_aq(
    R = Modelica.Constants.R/Mg_For2_aq.MM,
    G_ref = -1168872,
    H_ref = -1344215,
    S_ref = -2.966,
    a1 = 3.4516e-05,
    a2 = 5171.97,
    a3 = 3.7334e-05,
    a4 = -137649,
    c1 = 202.0613,
    c2 = 489687,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.050346) "Mg(For)2,aq";

  //1192
  constant DataRecord Ca_For_p1(
    R = Modelica.Constants.R/Ca_For_p1.MM,
    G_ref = -911832,
    H_ref = -970680,
    S_ref = 54.413,
    a1 = 1.6462e-05,
    a2 = 764.17,
    a3 = 0.00021048,
    a4 = -119428,
    c1 = 114.7847,
    c2 = 134901,
    w_ref = 148070,
    z = 1,
    a_0 = 3.72,
    MM = 0.085098) "Ca(For)+";

  //1193
  constant DataRecord Ca_For2_aq(
    R = Modelica.Constants.R/Ca_For2_aq.MM,
    G_ref = -1267677,
    H_ref = -1401849,
    S_ref = 142.607,
    a1 = 3.6997e-05,
    a2 = 5777.85,
    a3 = 1.3472e-05,
    a4 = -140156,
    c1 = 178.4869,
    c2 = 407748,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.066116) "Ca(For)2,aq";

  //1194
  constant DataRecord Sr_For_p1(
    R = Modelica.Constants.R/Sr_For_p1.MM,
    G_ref = -922647,
    H_ref = -975571,
    S_ref = 88.471,
    a1 = 1.6702e-05,
    a2 = 821.57,
    a3 = 0.00020845,
    a4 = -119667,
    c1 = 96.1027,
    c2 = 86575,
    w_ref = 96190,
    z = 1,
    a_0 = 3.72,
    MM = 0.13264) "Sr(For)+";

  //1195
  constant DataRecord Sr_For2_aq(
    R = Modelica.Constants.R/Sr_For2_aq.MM,
    G_ref = -1278266,
    H_ref = -1403376,
    S_ref = 187.175,
    a1 = 3.7459e-05,
    a2 = 5891.78,
    a3 = 8.749e-06,
    a4 = -140624,
    c1 = 151.36,
    c2 = 313461,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.11366) "Sr(For)2,aq";

  //1196
  constant DataRecord Ba_For_p1(
    R = Modelica.Constants.R/Ba_For_p1.MM,
    G_ref = -919539,
    H_ref = -957793,
    S_ref = 144.545,
    a1 = 1.9484e-05,
    a2 = 1502.31,
    a3 = 0.00018141,
    a4 = -122478,
    c1 = 75.9024,
    c2 = 43430,
    w_ref = 11670,
    z = 1,
    a_0 = 3.72,
    MM = 0.18236) "Ba(For)+";

  //1197
  constant DataRecord Ba_For2_aq(
    R = Modelica.Constants.R/Ba_For2_aq.MM,
    G_ref = -1275153,
    H_ref = -1380440,
    S_ref = 260.559,
    a1 = 4.0878e-05,
    a2 = 6724.48,
    a3 = -2.3527e-05,
    a4 = -144068,
    c1 = 127.1392,
    c2 = 229279,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16338) "Ba(For)2,aq";

  //1198
  constant DataRecord Cu_For_p1(
    R = Modelica.Constants.R/Cu_For_p1.MM,
    G_ref = -296595,
    H_ref = -369631,
    S_ref = -0.92,
    a1 = 1.2573e-05,
    a2 = -185.35,
    a3 = 0.00024774,
    a4 = -115504,
    c1 = 132.578,
    c2 = 169992,
    w_ref = 231630,
    z = 1,
    a_0 = 3.72,
    MM = 0.10856) "Cu(For)+";

  //1199
  constant DataRecord Cu_For2_aq(
    R = Modelica.Constants.R/Cu_For2_aq.MM,
    G_ref = -655009,
    H_ref = -808278,
    S_ref = 70.195,
    a1 = 3.2349e-05,
    a2 = 4643.36,
    a3 = 5.7953e-05,
    a4 = -135465,
    c1 = 198.186,
    c2 = 476219,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.089576) "Cu(For)2,aq";

  //1200
  constant DataRecord Cd_For_p1(
    R = Modelica.Constants.R/Cd_For_p1.MM,
    G_ref = -439035,
    H_ref = -507603,
    S_ref = 32.167,
    a1 = 1.8151e-05,
    a2 = 1175.58,
    a3 = 0.00019447,
    a4 = -121131,
    c1 = 140.2577,
    c2 = 212179,
    w_ref = 183220,
    z = 1,
    a_0 = 3.72,
    MM = 0.15742) "Cd(For)+";

  //1201
  constant DataRecord Cd_For2_aq(
    R = Modelica.Constants.R/Cd_For2_aq.MM,
    G_ref = -801332,
    H_ref = -947270,
    S_ref = 113.491,
    a1 = 3.8745e-05,
    a2 = 6204.83,
    a3 = -3.314e-06,
    a4 = -141921,
    c1 = 221.8683,
    c2 = 558531,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13844) "Cd(For)2,aq";

  //1202
  constant DataRecord Na_For_aq(
    R = Modelica.Constants.R/Na_For_aq.MM,
    G_ref = -613044,
    H_ref = -666423,
    S_ref = 147.833,
    a1 = 2.6732e-05,
    a2 = 3271.72,
    a3 = 0.00011195,
    a4 = -129796,
    c1 = 80.437,
    c2 = 66952,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.036008) "Na(For),aq";

  //1203
  constant DataRecord Na_For2_n1(
    R = Modelica.Constants.R/Na_For2_n1.MM,
    G_ref = -962324,
    H_ref = -1103425,
    S_ref = 190.259,
    a1 = 5.0396e-05,
    a2 = 9050.33,
    a3 = -0.00011534,
    a4 = -153683,
    c1 = 171.7524,
    c2 = 256299,
    w_ref = 387310,
    z = -1,
    a_0 = 3.72,
    MM = 0.11303) "Na(For)2-";

  //1204
  constant DataRecord K_For_aq(
    R = Modelica.Constants.R/K_For_aq.MM,
    G_ref = -633512,
    H_ref = -674256,
    S_ref = 203.627,
    a1 = 3.322e-05,
    a2 = 4856.37,
    a3 = 4.9534e-05,
    a4 = -136344,
    c1 = 41.3752,
    c2 = -68814,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.052118) "K(For),aq";

  //1205
  constant DataRecord K_For2_n1(
    R = Modelica.Constants.R/K_For2_n1.MM,
    G_ref = -982679,
    H_ref = -1106923,
    S_ref = 264.4,
    a1 = 5.7266e-05,
    a2 = 10728.11,
    a3 = -0.00018136,
    a4 = -160620,
    c1 = 85.7607,
    c2 = -8598,
    w_ref = 281210,
    z = -1,
    a_0 = 3.72,
    MM = 0.12914) "K(For)2-";

  //1206
  constant DataRecord La_For_p2(
    R = Modelica.Constants.R/La_For_p2.MM,
    G_ref = -1052067,
    H_ref = -1147851,
    S_ref = -119.403,
    a1 = 4.9614e-06,
    a2 = -2042.13,
    a3 = 0.00032036,
    a4 = -107826,
    c1 = 105.8339,
    c2 = -47974,
    w_ref = 622040,
    z = 2,
    a_0 = 3.72,
    MM = 0.18393) "La(For)+2";

  //1207
  constant DataRecord La_For2_p1(
    R = Modelica.Constants.R/La_For2_p1.MM,
    G_ref = -1413560,
    H_ref = -1584719,
    S_ref = -31.372,
    a1 = 2.338e-05,
    a2 = 2452.95,
    a3 = 0.00014416,
    a4 = -126411,
    c1 = 108.025,
    c2 = 69806,
    w_ref = 277980,
    z = 1,
    a_0 = 3.72,
    MM = 0.22895) "La(For)2+";

  //1208
  constant DataRecord Eu_For_p1(
    R = Modelica.Constants.R/Eu_For_p1.MM,
    G_ref = -899024,
    H_ref = -949994,
    S_ref = 117.734,
    a1 = 2.6257e-05,
    a2 = 3156.03,
    a3 = 0.00011637,
    a4 = -129315,
    c1 = 180.6241,
    c2 = 394355,
    w_ref = 52470,
    z = 1,
    a_0 = 3.72,
    MM = 0.19698) "Eu(For)+";

  //1209
  constant DataRecord Eu_For2_aq(
    R = Modelica.Constants.R/Eu_For2_aq.MM,
    G_ref = -1257380,
    H_ref = -1377850,
    S_ref = 225.472,
    a1 = 4.827e-05,
    a2 = 8529.17,
    a3 = -9.4412e-05,
    a4 = -151528,
    c1 = 324.1324,
    c2 = 913974,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.178) "Eu(For)2,aq";

  //1210
  constant DataRecord U_For_p2(
    R = Modelica.Constants.R/U_For_p2.MM,
    G_ref = -842695,
    H_ref = -925986,
    S_ref = -84.391,
    a1 = 4.8518e-06,
    a2 = -2068.23,
    a3 = 0.00032126,
    a4 = -107717,
    c1 = 102.7653,
    c2 = -42221,
    w_ref = 570780,
    z = 2,
    a_0 = 3.72,
    MM = 0.28305) "U(For)+2";

  //1211
  constant DataRecord U_For2_p1(
    R = Modelica.Constants.R/U_For2_p1.MM,
    G_ref = -1208302,
    H_ref = -1363653,
    S_ref = 16.364,
    a1 = 2.3208e-05,
    a2 = 2411.82,
    a3 = 0.00014559,
    a4 = -126240,
    c1 = 104.6293,
    c2 = 81032,
    w_ref = 206060,
    z = 1,
    a_0 = 3.72,
    MM = 0.32807) "U(For)2+";

  //1212
  constant DataRecord Eu_For_p2(
    R = Modelica.Constants.R/Eu_For_p2.MM,
    G_ref = -941266,
    H_ref = -1045105,
    S_ref = -125.336,
    a1 = 3.2673e-06,
    a2 = -2454.63,
    a3 = 0.00033633,
    a4 = -106123,
    c1 = 109.9137,
    c2 = -36468,
    w_ref = 630400,
    z = 2,
    a_0 = 3.72,
    MM = 0.19698) "Eu(For)+2";

  //1213
  constant DataRecord Eu_For2_p1(
    R = Modelica.Constants.R/Eu_For2_p1.MM,
    G_ref = -1303563,
    H_ref = -1483412,
    S_ref = -39.459,
    a1 = 2.1499e-05,
    a2 = 1992.84,
    a3 = 0.00016239,
    a4 = -124507,
    c1 = 115.514,
    c2 = 92257,
    w_ref = 289160,
    z = 1,
    a_0 = 3.72,
    MM = 0.242) "Eu(For)2+";

  //1214
  constant DataRecord Gd_For_p2(
    R = Modelica.Constants.R/Gd_For_p2.MM,
    G_ref = -1029645,
    H_ref = -1123965,
    S_ref = -102.788,
    a1 = 3.7309e-06,
    a2 = -2341.12,
    a3 = 0.00033182,
    a4 = -106592,
    c1 = 110.7689,
    c2 = -23045,
    w_ref = 597770,
    z = 2,
    a_0 = 3.72,
    MM = 0.20227) "Gd(For)+2";

  //1215
  constant DataRecord Gd_For2_p1(
    R = Modelica.Constants.R/Gd_For2_p1.MM,
    G_ref = -1391310,
    H_ref = -1559205,
    S_ref = -8.719,
    a1 = 2.1985e-05,
    a2 = 2113.38,
    a3 = 0.00015732,
    a4 = -125005,
    c1 = 118.8858,
    c2 = 118445,
    w_ref = 243970,
    z = 1,
    a_0 = 3.72,
    MM = 0.24729) "Gd(For)2+";

  //1216
  constant DataRecord Yb_For_p2(
    R = Modelica.Constants.R/Yb_For_p2.MM,
    G_ref = -1006440,
    H_ref = -1111894,
    S_ref = -148.478,
    a1 = 1.3585e-06,
    a2 = -2920.56,
    a3 = 0.0003546,
    a4 = -104194,
    c1 = 114.6316,
    c2 = -32631,
    w_ref = 669650,
    z = 2,
    a_0 = 3.72,
    MM = 0.21806) "Yb(For)+2";

  //1217
  constant DataRecord Yb_For2_p1(
    R = Modelica.Constants.R/Yb_For2_p1.MM,
    G_ref = -1368281,
    H_ref = -1552256,
    S_ref = -71.015,
    a1 = 1.9394e-05,
    a2 = 1479.13,
    a3 = 0.00018254,
    a4 = -122382,
    c1 = 122.292,
    c2 = 99738,
    w_ref = 339360,
    z = 1,
    a_0 = 3.72,
    MM = 0.26308) "Yb(For)2+";

  //1218
  constant DataRecord Pb_For_p1(
    R = Modelica.Constants.R/Pb_For_p1.MM,
    G_ref = -385556,
    H_ref = -421279,
    S_ref = 155.385,
    a1 = 1.7517e-05,
    a2 = 1020.44,
    a3 = 0.00020067,
    a4 = -120487,
    c1 = 72.2355,
    c2 = 35756,
    w_ref = -4180,
    z = 1,
    a_0 = 3.72,
    MM = 0.25221) "Pb(For)+";

  //1219
  constant DataRecord Pb_For2_aq(
    R = Modelica.Constants.R/Pb_For2_aq.MM,
    G_ref = -743572,
    H_ref = -845327,
    S_ref = 274.742,
    a1 = 3.8745e-05,
    a2 = 6204.83,
    a3 = -3.314e-06,
    a4 = -141921,
    c1 = 122.8335,
    c2 = 214313,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.23323) "Pb(For)2,aq";

  //1220
  constant DataRecord Mn_For_p1(
    R = Modelica.Constants.R/Mn_For_p1.MM,
    G_ref = -588609,
    H_ref = -651595,
    S_ref = 31.024,
    a1 = 1.7194e-05,
    a2 = 943.28,
    a3 = 0.0002033,
    a4 = -120169,
    c1 = 136.9473,
    c2 = 200673,
    w_ref = 183220,
    z = 1,
    a_0 = 3.72,
    MM = 0.099958) "Mn(For)+";

  //1221
  constant DataRecord Mn_For2_aq(
    R = Modelica.Constants.R/Mn_For2_aq.MM,
    G_ref = -945710,
    H_ref = -1086171,
    S_ref = 112.001,
    a1 = 3.7679e-05,
    a2 = 5943.79,
    a3 = 7.117e-06,
    a4 = -140842,
    c1 = 215.4099,
    c2 = 536079,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.080976) "Mn(For)2,aq";

  //1222
  constant DataRecord Co_For_p1(
    R = Modelica.Constants.R/Co_For_p1.MM,
    G_ref = -415885,
    H_ref = -494331,
    S_ref = -22.598,
    a1 = 1.2809e-05,
    a2 = -127.07,
    a3 = 0.00024534,
    a4 = -115742,
    c1 = 123.9569,
    c2 = 129721,
    w_ref = 263800,
    z = 1,
    a_0 = 3.72,
    MM = 0.10395) "Co(For)+";

  //1223
  constant DataRecord Co_For2_aq(
    R = Modelica.Constants.R/Co_For2_aq.MM,
    G_ref = -773730,
    H_ref = -934584,
    S_ref = 41.827,
    a1 = 3.2491e-05,
    a2 = 4676.58,
    a3 = 5.6965e-05,
    a4 = -135603,
    c1 = 175.5811,
    c2 = 397643,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.084966) "Co(For)2,aq";

  //1224
  constant DataRecord Ni_For_p1(
    R = Modelica.Constants.R/Ni_For_p1.MM,
    G_ref = -407158,
    H_ref = -491925,
    S_ref = -44.275,
    a1 = 9.9864e-06,
    a2 = -818.31,
    a3 = 0.00027295,
    a4 = -112884,
    c1 = 105.4891,
    c2 = 54936,
    w_ref = 296900,
    z = 1,
    a_0 = 3.72,
    MM = 0.10373) "Ni(For)+";

  //1225
  constant DataRecord Ni_For2_aq(
    R = Modelica.Constants.R/Ni_For2_aq.MM,
    G_ref = -765057,
    H_ref = -934668,
    S_ref = 13.46,
    a1 = 2.9221e-05,
    a2 = 3878.02,
    a3 = 8.842e-05,
    a4 = -132302,
    c1 = 133.5985,
    c2 = 251726,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.084746) "Ni(For)2,aq";

  //1226
  constant DataRecord Zn_For_p1(
    R = Modelica.Constants.R/Zn_For_p1.MM,
    G_ref = -508260,
    H_ref = -588680,
    S_ref = -18.033,
    a1 = 1.285e-05,
    a2 = -116.94,
    a3 = 0.0002449,
    a4 = -115784,
    c1 = 131.6065,
    c2 = 158486,
    w_ref = 257020,
    z = 1,
    a_0 = 3.72,
    MM = 0.11039) "Zn(For)+";

  //1227
  constant DataRecord Zn_For2_aq(
    R = Modelica.Constants.R/Zn_For2_aq.MM,
    G_ref = -865703,
    H_ref = -1028118,
    S_ref = 47.798,
    a1 = 3.2562e-05,
    a2 = 4694.62,
    a3 = 5.6137e-05,
    a4 = -135679,
    c1 = 191.7276,
    c2 = 453767,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.091406) "Zn(For)2,aq";

  //1228
  constant DataRecord Fe_For_p1(
    R = Modelica.Constants.R/Fe_For_p1.MM,
    G_ref = -452813,
    H_ref = -525724,
    S_ref = -14.611,
    a1 = 1.4178e-05,
    a2 = 205.98,
    a3 = 0.00023253,
    a4 = -117123,
    c1 = 121.3695,
    c2 = 123968,
    w_ref = 253680,
    z = 1,
    a_0 = 3.72,
    MM = 0.10087) "Fe(For)+";

  //1229
  constant DataRecord Fe_For2_aq(
    R = Modelica.Constants.R/Fe_For2_aq.MM,
    G_ref = -810483,
    H_ref = -965073,
    S_ref = 52.279,
    a1 = 3.4054e-05,
    a2 = 5060.63,
    a3 = 4.1388e-05,
    a4 = -137189,
    c1 = 172.3511,
    c2 = 386422,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.081886) "Fe(For)2,aq";

  //1230
  constant DataRecord Mg_Prop_p1(
    R = Modelica.Constants.R/Mg_Prop_p1.MM,
    G_ref = -821139,
    H_ref = -986001,
    S_ref = -36.744,
    a1 = 3.1183e-05,
    a2 = 4358.6,
    a3 = 6.9233e-05,
    a4 = -134290,
    c1 = 380.2683,
    c2 = 1013687,
    w_ref = 285390,
    z = 1,
    a_0 = 3.72,
    MM = 0.097379) "Mg(Prop)+";

  //1231
  constant DataRecord Mg_Prop2_aq(
    R = Modelica.Constants.R/Mg_Prop2_aq.MM,
    G_ref = -1185900,
    H_ref = -1509960,
    S_ref = 43.401,
    a1 = 6.9365e-05,
    a2 = 13680.42,
    a3 = -0.00029692,
    a4 = -172824,
    c1 = 776.9592,
    c2 = 2487882,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.10645) "Mg(Prop)2,aq";

  //1232
  constant DataRecord Ca_Prop_p1(
    R = Modelica.Constants.R/Ca_Prop_p1.MM,
    G_ref = -919719,
    H_ref = -1054054,
    S_ref = 74.496,
    a1 = 3.2846e-05,
    a2 = 4763.78,
    a3 = 5.3455e-05,
    a4 = -135963,
    c1 = 352.8112,
    c2 = 971688,
    w_ref = 118490,
    z = 1,
    a_0 = 3.72,
    MM = 0.11315) "Ca(Prop)+";

  //1233
  constant DataRecord Ca_Prop2_aq(
    R = Modelica.Constants.R/Ca_Prop2_aq.MM,
    G_ref = -1284245,
    H_ref = -1567548,
    S_ref = 188.97,
    a1 = 7.1846e-05,
    a2 = 14286.31,
    a3 = -0.00032078,
    a4 = -175331,
    c1 = 753.3848,
    c2 = 2405942,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.12222) "Ca(Prop)2,aq";

  //1234
  constant DataRecord Sr_Prop_p1(
    R = Modelica.Constants.R/Sr_Prop_p1.MM,
    G_ref = -928250,
    H_ref = -1056661,
    S_ref = 108.554,
    a1 = 3.3085e-05,
    a2 = 4823.9,
    a3 = 5.0693e-05,
    a4 = -136210,
    c1 = 334.115,
    c2 = 923367,
    w_ref = 66480,
    z = 1,
    a_0 = 3.72,
    MM = 0.16069) "Sr(Prop)+";

  //1235
  constant DataRecord Sr_Prop2_aq(
    R = Modelica.Constants.R/Sr_Prop2_aq.MM,
    G_ref = -1290726,
    H_ref = -1564967,
    S_ref = 233.538,
    a1 = 7.2308e-05,
    a2 = 14400.24,
    a3 = -0.0003255,
    a4 = -175799,
    c1 = 726.2579,
    c2 = 2311656,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16976) "Sr(Prop)2,aq";

  //1236
  constant DataRecord Ba_Prop_p1(
    R = Modelica.Constants.R/Ba_Prop_p1.MM,
    G_ref = -924685,
    H_ref = -1038427,
    S_ref = 164.628,
    a1 = 3.5865e-05,
    a2 = 5502.59,
    a3 = 2.4083e-05,
    a4 = -139018,
    c1 = 313.859,
    c2 = 880217,
    w_ref = -18660,
    z = 1,
    a_0 = 3.72,
    MM = 0.21041) "Ba(Prop)+";

  //1237
  constant DataRecord Ba_Prop2_aq(
    R = Modelica.Constants.R/Ba_Prop2_aq.MM,
    G_ref = -1286701,
    H_ref = -1541118,
    S_ref = 306.922,
    a1 = 7.5727e-05,
    a2 = 15235.53,
    a3 = -0.00035845,
    a4 = -179251,
    c1 = 702.0371,
    c2 = 2227474,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.21948) "Ba(Prop)2,aq";

  //1238
  constant DataRecord Eu_Prop_p1(
    R = Modelica.Constants.R/Eu_Prop_p1.MM,
    G_ref = -906455,
    H_ref = -1032912,
    S_ref = 137.817,
    a1 = 4.2636e-05,
    a2 = 7154.05,
    a3 = -4.0413e-05,
    a4 = -145846,
    c1 = 418.551,
    c2 = 1231146,
    w_ref = 21800,
    z = 1,
    a_0 = 3.72,
    MM = 0.22503) "Eu(Prop)+";

  //1239
  constant DataRecord Eu_Prop2_aq(
    R = Modelica.Constants.R/Eu_Prop2_aq.MM,
    G_ref = -1268530,
    H_ref = -1538126,
    S_ref = 271.834,
    a1 = 8.3118e-05,
    a2 = 17040.18,
    a3 = -0.0004293,
    a4 = -186715,
    c1 = 899.0303,
    c2 = 2912169,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.2341) "Eu(Prop)2,aq";

  //1240
  constant DataRecord Cu_Prop_p1(
    R = Modelica.Constants.R/Cu_Prop_p1.MM,
    G_ref = -310135,
    H_ref = -458470,
    S_ref = 19.163,
    a1 = 2.8952e-05,
    a2 = 3812.71,
    a3 = 9.0872e-05,
    a4 = -132030,
    c1 = 370.4873,
    c2 = 1006783,
    w_ref = 200790,
    z = 1,
    a_0 = 3.72,
    MM = 0.13661) "Cu(Prop)+";

  //1241
  constant DataRecord Cu_Prop2_aq(
    R = Modelica.Constants.R/Cu_Prop2_aq.MM,
    G_ref = -681971,
    H_ref = -984361,
    S_ref = 116.558,
    a1 = 6.7198e-05,
    a2 = 13151.82,
    a3 = -0.0002763,
    a4 = -170640,
    c1 = 773.0835,
    c2 = 2474413,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14568) "Cu(Prop)2,aq";

  //1242
  constant DataRecord Cd_Prop_p1(
    R = Modelica.Constants.R/Cd_Prop_p1.MM,
    G_ref = -451487,
    H_ref = -595542,
    S_ref = 52.25,
    a1 = 3.4529e-05,
    a2 = 5173.93,
    a3 = 3.7476e-05,
    a4 = -137658,
    c1 = 378.1441,
    c2 = 1048971,
    w_ref = 152130,
    z = 1,
    a_0 = 3.72,
    MM = 0.18547) "Cd(Prop)+";

  //1243
  constant DataRecord Cd_Prop2_aq(
    R = Modelica.Constants.R/Cd_Prop2_aq.MM,
    G_ref = -822240,
    H_ref = -1117308,
    S_ref = 159.858,
    a1 = 7.3594e-05,
    a2 = 14713.25,
    a3 = -0.00033757,
    a4 = -177092,
    c1 = 796.7662,
    c2 = 2556725,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19454) "Cd(Prop)2,aq";

  //1244
  constant DataRecord Pb_Prop_p1(
    R = Modelica.Constants.R/Pb_Prop_p1.MM,
    G_ref = -400292,
    H_ref = -511502,
    S_ref = 175.469,
    a1 = 3.3896e-05,
    a2 = 5021.26,
    a3 = 4.3108e-05,
    a4 = -137026,
    c1 = 310.1386,
    c2 = 872548,
    w_ref = -35060,
    z = 1,
    a_0 = 3.72,
    MM = 0.28026) "Pb(Prop)+";

  //1245
  constant DataRecord Pb_Prop2_aq(
    R = Modelica.Constants.R/Pb_Prop2_aq.MM,
    G_ref = -770701,
    H_ref = -1021582,
    S_ref = 321.105,
    a1 = 7.3594e-05,
    a2 = 14713.25,
    a3 = -0.00033757,
    a4 = -177092,
    c1 = 697.7314,
    c2 = 2212508,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.28933) "Pb(Prop)2,aq";

  //1246
  constant DataRecord Na_Prop_aq(
    R = Modelica.Constants.R/Na_Prop_aq.MM,
    G_ref = -625211,
    H_ref = -754078,
    S_ref = 167.916,
    a1 = 4.3216e-05,
    a2 = 7297.23,
    a3 = -4.6442e-05,
    a4 = -146436,
    c1 = 321.1889,
    c2 = 903740,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.06406) "Na(Prop),aq";

  //1247
  constant DataRecord Na_Prop2_n1(
    R = Modelica.Constants.R/Na_Prop2_n1.MM,
    G_ref = -986775,
    H_ref = -1277329,
    S_ref = 239.706,
    a1 = 8.5013e-05,
    a2 = 17502.13,
    a3 = -0.00044737,
    a4 = -188623,
    c1 = 740.3224,
    c2 = 2254494,
    w_ref = 318650,
    z = -1,
    a_0 = 3.72,
    MM = 0.16913) "Na(Prop)2-";

  //1248
  constant DataRecord K_Prop_aq(
    R = Modelica.Constants.R/K_Prop_aq.MM,
    G_ref = -645679,
    H_ref = -761911,
    S_ref = 223.71,
    a1 = 4.9703e-05,
    a2 = 8879.41,
    a3 = -0.0001082,
    a4 = -152975,
    c1 = 282.1267,
    c2 = 767973,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.08017) "K(Prop),aq";

  //1249
  constant DataRecord K_Prop2_n1(
    R = Modelica.Constants.R/K_Prop2_n1.MM,
    G_ref = -1007126,
    H_ref = -1280827,
    S_ref = 309.662,
    a1 = 9.1883e-05,
    a2 = 19179.79,
    a3 = -0.00051332,
    a4 = -195560,
    c1 = 654.3441,
    c2 = 1989597,
    w_ref = 212630,
    z = -1,
    a_0 = 3.72,
    MM = 0.18524) "K(Prop)2-";

  //1250
  constant DataRecord La_Prop_p2(
    R = Modelica.Constants.R/La_Prop_p2.MM,
    G_ref = -1063606,
    H_ref = -1234874,
    S_ref = -99.32,
    a1 = 2.1349e-05,
    a2 = 1956.27,
    a3 = 0.00016388,
    a4 = -124357,
    c1 = 343.9855,
    c2 = 788810,
    w_ref = 593790,
    z = 2,
    a_0 = 3.72,
    MM = 0.21198) "La(Prop)+2";

  //1251
  constant DataRecord La_Prop2_p1(
    R = Modelica.Constants.R/La_Prop2_p1.MM,
    G_ref = -1436869,
    H_ref = -1756824,
    S_ref = 16.096,
    a1 = 5.7986e-05,
    a2 = 10902.37,
    a3 = -0.00018785,
    a4 = -161339,
    c1 = 676.2967,
    c2 = 2068005,
    w_ref = 206060,
    z = 1,
    a_0 = 3.72,
    MM = 0.28505) "La(Prop)2+";

  //1252
  constant DataRecord Eu_Prop_p2(
    R = Modelica.Constants.R/Eu_Prop_p2.MM,
    G_ref = -953835,
    H_ref = -1133157,
    S_ref = -105.253,
    a1 = 1.9654e-05,
    a2 = 1544.11,
    a3 = 0.00017967,
    a4 = -122654,
    c1 = 348.0251,
    c2 = 800316,
    w_ref = 601740,
    z = 2,
    a_0 = 3.72,
    MM = 0.22503) "Eu(Prop)+2";

  //1253
  constant DataRecord Eu_Prop2_p1(
    R = Modelica.Constants.R/Eu_Prop2_p1.MM,
    G_ref = -1328696,
    H_ref = -1657345,
    S_ref = 8.004,
    a1 = 5.6113e-05,
    a2 = 10445.19,
    a3 = -0.00016988,
    a4 = -159448,
    c1 = 684.03,
    c2 = 2090452,
    w_ref = 219910,
    z = 1,
    a_0 = 3.72,
    MM = 0.2981) "Eu(Prop)2+";

  //1254
  constant DataRecord Gd_Prop_p2(
    R = Modelica.Constants.R/Gd_Prop_p2.MM,
    G_ref = -1042155,
    H_ref = -1211963,
    S_ref = -82.705,
    a1 = 2.011e-05,
    a2 = 1654.23,
    a3 = 0.00017566,
    a4 = -123106,
    c1 = 348.6895,
    c2 = 813742,
    w_ref = 567020,
    z = 2,
    a_0 = 3.72,
    MM = 0.23032) "Gd(Prop)+2";

  //1255
  constant DataRecord Gd_Prop2_p1(
    R = Modelica.Constants.R/Gd_Prop2_p1.MM,
    G_ref = -1416389,
    H_ref = -1733080,
    S_ref = 38.748,
    a1 = 5.659e-05,
    a2 = 10563.14,
    a3 = -0.00017487,
    a4 = -159938,
    c1 = 687.1111,
    c2 = 2116644,
    w_ref = 171540,
    z = 1,
    a_0 = 3.72,
    MM = 0.30339) "Gd(Prop)2+";

  //1256
  constant DataRecord Yb_Prop_p2(
    R = Modelica.Constants.R/Yb_Prop_p2.MM,
    G_ref = -1017867,
    H_ref = -1198808,
    S_ref = -128.394,
    a1 = 1.7738e-05,
    a2 = 1074.83,
    a3 = 0.00019843,
    a4 = -120713,
    c1 = 352.5505,
    c2 = 804152,
    w_ref = 638850,
    z = 2,
    a_0 = 3.72,
    MM = 0.24611) "Yb(Prop)+2";

  //1257
  constant DataRecord Yb_Prop2_p1(
    R = Modelica.Constants.R/Yb_Prop2_p1.MM,
    G_ref = -1391360,
    H_ref = -1724134,
    S_ref = -23.548,
    a1 = 5.3999e-05,
    a2 = 9928.63,
    a3 = -0.00014951,
    a4 = -157314,
    c1 = 690.547,
    c2 = 2097937,
    w_ref = 267270,
    z = 1,
    a_0 = 3.72,
    MM = 0.31918) "Yb(Prop)2+";

  //1258
  constant DataRecord U_Prop_p2(
    R = Modelica.Constants.R/U_Prop_p2.MM,
    G_ref = -854862,
    H_ref = -1013641,
    S_ref = -64.308,
    a1 = 2.1237e-05,
    a2 = 1930.92,
    a3 = 0.0001644,
    a4 = -124252,
    c1 = 340.8291,
    c2 = 794563,
    w_ref = 541580,
    z = 2,
    a_0 = 3.72,
    MM = 0.3111) "U(Prop)+2";

  //1259
  constant DataRecord U_Prop2_p1(
    R = Modelica.Constants.R/U_Prop2_p1.MM,
    G_ref = -1228125,
    H_ref = -1531796,
    S_ref = 63.827,
    a1 = 5.7815e-05,
    a2 = 10860.91,
    a3 = -0.00018626,
    a4 = -161168,
    c1 = 672.9362,
    c2 = 2079226,
    w_ref = 134520,
    z = 1,
    a_0 = 3.72,
    MM = 0.38417) "U(Prop)2+";

  //1260
  constant DataRecord Co_Prop_p1(
    R = Modelica.Constants.R/Co_Prop_p1.MM,
    G_ref = -424914,
    H_ref = -578844,
    S_ref = -2.515,
    a1 = 2.9194e-05,
    a2 = 3872.08,
    a3 = 8.85e-05,
    a4 = -132277,
    c1 = 362.0214,
    c2 = 966512,
    w_ref = 234640,
    z = 1,
    a_0 = 3.72,
    MM = 0.132) "Co(Prop)+";

  //1261
  constant DataRecord Co_Prop2_aq(
    R = Modelica.Constants.R/Co_Prop2_aq.MM,
    G_ref = -792471,
    H_ref = -1102451,
    S_ref = 88.19,
    a1 = 6.734e-05,
    a2 = 13187.63,
    a3 = -0.00027796,
    a4 = -170787,
    c1 = 750.4778,
    c2 = 2395842,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14107) "Co(Prop)2,aq";

  //1262
  constant DataRecord Ni_Prop_p1(
    R = Modelica.Constants.R/Ni_Prop_p1.MM,
    G_ref = -416873,
    H_ref = -577124,
    S_ref = -24.192,
    a1 = 2.637e-05,
    a2 = 3183.9,
    a3 = 0.00011525,
    a4 = -129432,
    c1 = 343.5106,
    c2 = 891723,
    w_ref = 267270,
    z = 1,
    a_0 = 3.72,
    MM = 0.13178) "Ni(Prop)+";

  //1263
  constant DataRecord Ni_Prop2_aq(
    R = Modelica.Constants.R/Ni_Prop2_aq.MM,
    G_ref = -785052,
    H_ref = -1103354,
    S_ref = 59.823,
    a1 = 6.407e-05,
    a2 = 12389.03,
    a3 = -0.00024648,
    a4 = -167486,
    c1 = 708.4964,
    c2 = 2249921,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14085) "Ni(Prop)2,aq";

  //1264
  constant DataRecord Zn_Prop_p1(
    R = Modelica.Constants.R/Zn_Prop_p1.MM,
    G_ref = -517456,
    H_ref = -673369,
    S_ref = 2.05,
    a1 = 2.9237e-05,
    a2 = 3884.17,
    a3 = 8.7705e-05,
    a4 = -132327,
    c1 = 369.7434,
    c2 = 995277,
    w_ref = 228610,
    z = 1,
    a_0 = 3.72,
    MM = 0.13844) "Zn(Prop)+";

  //1265
  constant DataRecord Zn_Prop2_aq(
    R = Modelica.Constants.R/Zn_Prop2_aq.MM,
    G_ref = -884728,
    H_ref = -1196268,
    S_ref = 94.165,
    a1 = 6.7411e-05,
    a2 = 13202.95,
    a3 = -0.00027812,
    a4 = -170849,
    c1 = 766.6255,
    c2 = 2451962,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.14751) "Zn(Prop)2,aq";

  //1266
  constant DataRecord Fe_Prop_p1(
    R = Modelica.Constants.R/Fe_Prop_p1.MM,
    G_ref = -463725,
    H_ref = -614450,
    S_ref = 5.473,
    a1 = 3.0557e-05,
    a2 = 4206.68,
    a3 = 7.4952e-05,
    a4 = -133658,
    c1 = 359.2725,
    c2 = 960759,
    w_ref = 222760,
    z = 1,
    a_0 = 3.72,
    MM = 0.12892) "Fe(Prop)+";

  //1267
  constant DataRecord Fe_Prop2_aq(
    R = Modelica.Constants.R/Fe_Prop2_aq.MM,
    G_ref = -832649,
    H_ref = -1138784,
    S_ref = 98.642,
    a1 = 6.8903e-05,
    a2 = 13569.09,
    a3 = -0.00029286,
    a4 = -172364,
    c1 = 747.2486,
    c2 = 2384616,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13799) "Fe(Prop)2,aq";

  //1268
  constant DataRecord Mn_Prop_p1(
    R = Modelica.Constants.R/Mn_Prop_p1.MM,
    G_ref = -598379,
    H_ref = -736853,
    S_ref = 51.108,
    a1 = 3.3579e-05,
    a2 = 4942.48,
    a3 = 4.6488e-05,
    a4 = -136700,
    c1 = 375.0241,
    c2 = 1037465,
    w_ref = 154220,
    z = 1,
    a_0 = 3.72,
    MM = 0.12801) "Mn(Prop)+";

  //1269
  constant DataRecord Mn_Prop2_aq(
    R = Modelica.Constants.R/Mn_Prop2_aq.MM,
    G_ref = -965763,
    H_ref = -1255355,
    S_ref = 158.364,
    a1 = 7.2528e-05,
    a2 = 14454.72,
    a3 = -0.0003278,
    a4 = -176025,
    c1 = 790.307,
    c2 = 2534278,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.13708) "Mn(Prop)2,aq";

  //1270
  constant DataRecord U_But_p2(
    R = Modelica.Constants.R/U_But_p2.MM,
    G_ref = -845825,
    H_ref = -1038795,
    S_ref = -42.133,
    a1 = 2.9909e-05,
    a2 = 4047.52,
    a3 = 8.145e-05,
    a4 = -133001,
    c1 = 399.5699,
    c2 = 1009695,
    w_ref = 507350,
    z = 2,
    a_0 = 3.72,
    MM = 0.32513) "U(But)+2";

  //1271
  constant DataRecord U_But2_p1(
    R = Modelica.Constants.R/U_But2_p1.MM,
    G_ref = -1211423,
    H_ref = -1581037,
    S_ref = 116.24,
    a1 = 7.6126e-05,
    a2 = 15331.98,
    a3 = -0.00036201,
    a4 = -179653,
    c1 = 813.3817,
    c2 = 2592942,
    w_ref = 54730,
    z = 1,
    a_0 = 3.72,
    MM = 0.41222) "U(But)2+";

  //1272
  constant DataRecord Eu_But_p1(
    R = Modelica.Constants.R/Eu_But_p1.MM,
    G_ref = -895874,
    H_ref = -1057159,
    S_ref = 159.992,
    a1 = 5.1312e-05,
    a2 = 9272.5,
    a3 = -0.00012372,
    a4 = -154603,
    c1 = 477.3622,
    c2 = 1446275,
    w_ref = -11630,
    z = 1,
    a_0 = 3.72,
    MM = 0.23906) "Eu(But)+";

  //1273
  constant DataRecord Eu_But2_aq(
    R = Modelica.Constants.R/Eu_But2_aq.MM,
    G_ref = -1249029,
    H_ref = -1585560,
    S_ref = 323.03,
    a1 = 0.0001017,
    a2 = 21576.3,
    a3 = -0.00060747,
    a4 = -205464,
    c1 = 1046.8297,
    c2 = 3425884,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.26215) "Eu(But)2,aq";

  //1274
  constant DataRecord Mg_But_p1(
    R = Modelica.Constants.R/Mg_But_p1.MM,
    G_ref = -811190,
    H_ref = -1010177,
    S_ref = -14.569,
    a1 = 3.9865e-05,
    a2 = 6478.13,
    a3 = -1.4033e-05,
    a4 = -143051,
    c1 = 439.2422,
    c2 = 1228816,
    w_ref = 253680,
    z = 1,
    a_0 = 3.72,
    MM = 0.11141) "Mg(But)+";

  //1275
  constant DataRecord Mg_But2_aq(
    R = Modelica.Constants.R/Mg_But2_aq.MM,
    G_ref = -1166056,
    H_ref = -15563593,
    S_ref = 94.558,
    a1 = 8.7946e-05,
    a2 = 18219.14,
    a3 = -0.00047573,
    a4 = -191585,
    c1 = 924.7586,
    c2 = 3001597,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1345) "Mg(But)2,aq";

  //1276
  constant DataRecord Ca_But_p1(
    R = Modelica.Constants.R/Ca_But_p1.MM,
    G_ref = -909882,
    H_ref = -1078384,
    S_ref = 96.671,
    a1 = 4.152e-05,
    a2 = 6882.55,
    a3 = -2.9983e-05,
    a4 = -144720,
    c1 = 411.5931,
    c2 = 1186821,
    w_ref = 84730,
    z = 1,
    a_0 = 3.72,
    MM = 0.12718) "Ca(But)+";

  //1277
  constant DataRecord Ca_But2_aq(
    R = Modelica.Constants.R/Ca_But2_aq.MM,
    G_ref = -1264689,
    H_ref = -1614271,
    S_ref = 240.166,
    a1 = 9.0426e-05,
    a2 = 18825.11,
    a3 = -0.00049959,
    a4 = -194092,
    c1 = 901.1842,
    c2 = 2919658,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.15027) "Ca(But)2,aq";

  //1278
  constant DataRecord Sr_But_p1(
    R = Modelica.Constants.R/Sr_But_p1.MM,
    G_ref = -918702,
    H_ref = -1081141,
    S_ref = 130.729,
    a1 = 4.1759e-05,
    a2 = 6939.96,
    a3 = -3.2012e-05,
    a4 = -144959,
    c1 = 392.9107,
    c2 = 1138496,
    w_ref = 32840,
    z = 1,
    a_0 = 3.72,
    MM = 0.17472) "Sr(But)+";

  //1279
  constant DataRecord Sr_But2_aq(
    R = Modelica.Constants.R/Sr_But2_aq.MM,
    G_ref = -1271622,
    H_ref = -1612003,
    S_ref = 284.734,
    a1 = 9.0888e-05,
    a2 = 18936.45,
    a3 = -0.00050364,
    a4 = -194552,
    c1 = 874.0573,
    c2 = 2825372,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.20083) "Sr(But)2,aq";

  //1280
  constant DataRecord Ba_But_p1(
    R = Modelica.Constants.R/Ba_But_p1.MM,
    G_ref = -914790,
    H_ref = -1062610,
    S_ref = 186.803,
    a1 = 4.454e-05,
    a2 = 7621.16,
    a3 = -5.9245e-05,
    a4 = -147775,
    c1 = 372.6664,
    c2 = 1095346,
    w_ref = -52170,
    z = 1,
    a_0 = 3.72,
    MM = 0.22444) "Ba(But)+";

  //1281
  constant DataRecord Ba_But2_aq(
    R = Modelica.Constants.R/Ba_But2_aq.MM,
    G_ref = -1267028,
    H_ref = -1587619,
    S_ref = 358.117,
    a1 = 9.4307e-05,
    a2 = 19771.74,
    a3 = -0.00053659,
    a4 = -198004,
    c1 = 849.8373,
    c2 = 2741185,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.24753) "Ba(But)2,aq";

  //1282
  constant DataRecord Cu_But_p1(
    R = Modelica.Constants.R/Cu_But_p1.MM,
    G_ref = -300637,
    H_ref = -483093,
    S_ref = 41.338,
    a1 = 3.7634e-05,
    a2 = 5934.75,
    a3 = 7.037e-06,
    a4 = -140804,
    c1 = 429.4817,
    c2 = 1221912,
    w_ref = 169280,
    z = 1,
    a_0 = 3.72,
    MM = 0.15064) "Cu(But)+";

  //1283
  constant DataRecord Cu_But2_aq(
    R = Modelica.Constants.R/Cu_But2_aq.MM,
    G_ref = -663038,
    H_ref = -1031666,
    S_ref = 167.753,
    a1 = 8.5778e-05,
    a2 = 17688.03,
    a3 = -0.00045444,
    a4 = -189393,
    c1 = 920.8829,
    c2 = 2988129,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17373) "Cu(But)2,aq";

  //1284
  constant DataRecord Cd_But_p1(
    R = Modelica.Constants.R/Cd_But_p1.MM,
    G_ref = -440512,
    H_ref = -618663,
    S_ref = 74.425,
    a1 = 4.3204e-05,
    a2 = 7292.59,
    a3 = -4.5898e-05,
    a4 = -146415,
    c1 = 436.9397,
    c2 = 1264099,
    w_ref = 118490,
    z = 1,
    a_0 = 3.72,
    MM = 0.1995) "Cd(But)+";

  //1285
  constant DataRecord Cd_But2_aq(
    R = Modelica.Constants.R/Cd_But2_aq.MM,
    G_ref = -801081,
    H_ref = -1162353,
    S_ref = 211.054,
    a1 = 9.2175e-05,
    a2 = 19252.09,
    a3 = -0.00051638,
    a4 = -195857,
    c1 = 944.5656,
    c2 = 3070441,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22259) "Cd(But)2,aq";

  //1286
  constant DataRecord Na_But_aq(
    R = Modelica.Constants.R/Na_But_aq.MM,
    G_ref = -616173,
    H_ref = -776253,
    S_ref = 190.092,
    a1 = 5.2004e-05,
    a2 = 9440.94,
    a3 = -0.00013024,
    a4 = -155298,
    c1 = 383.0829,
    c2 = 1118873,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.078085) "Na(But),aq";

  //1287
  constant DataRecord Na_But2_n1(
    R = Modelica.Constants.R/Na_But2_n1.MM,
    G_ref = -968642,
    H_ref = -1325805,
    S_ref = 289.683,
    a1 = 0.00010334,
    a2 = 21977.38,
    a3 = -0.00062346,
    a4 = -207125,
    c1 = 881.145,
    c2 = 2768210,
    w_ref = 242920,
    z = -1,
    a_0 = 3.72,
    MM = 0.19718) "Na(But)2-";

  //1288
  constant DataRecord K_But_aq(
    R = Modelica.Constants.R/K_But_aq.MM,
    G_ref = -636642,
    H_ref = -779174,
    S_ref = 245.885,
    a1 = 5.8492e-05,
    a2 = 11025.68,
    a3 = -0.00019263,
    a4 = -161850,
    c1 = 344.0215,
    c2 = 983102,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.094195) "K(But),aq";

  //1289
  constant DataRecord K_But2_n1(
    R = Modelica.Constants.R/K_But2_n1.MM,
    G_ref = -988993,
    H_ref = -1329290,
    S_ref = 359.64,
    a1 = 0.00011021,
    a2 = 23655.04,
    a3 = -0.00068943,
    a4 = -214058,
    c1 = 795.1659,
    c2 = 2503308,
    w_ref = 136900,
    z = -1,
    a_0 = 3.72,
    MM = 0.21329) "K(But)2-";

  //1290
  constant DataRecord La_But_p2(
    R = Modelica.Constants.R/La_But_p2.MM,
    G_ref = -1055196,
    H_ref = -1260576,
    S_ref = -77.145,
    a1 = 3.0022e-05,
    a2 = 4075.38,
    a3 = 8.0283e-05,
    a4 = -133118,
    c1 = 402.7297,
    c2 = 1003942,
    w_ref = 559610,
    z = 2,
    a_0 = 3.72,
    MM = 0.22601) "La(But)+2";

  //1291
  constant DataRecord La_But2_p1(
    R = Modelica.Constants.R/La_But2_p1.MM,
    G_ref = -1419878,
    H_ref = -1805693,
    S_ref = 68.505,
    a1 = 7.63e-05,
    a2 = 15375.15,
    a3 = -0.00036387,
    a4 = -179828,
    c1 = 816.8331,
    c2 = 2581716,
    w_ref = 127240,
    z = 1,
    a_0 = 3.72,
    MM = 0.3131) "La(But)2+";

  //1292
  constant DataRecord Eu_But_p2(
    R = Modelica.Constants.R/Eu_But_p2.MM,
    G_ref = -944396,
    H_ref = -1158671,
    S_ref = -83.078,
    a1 = 2.8325e-05,
    a2 = 3661.13,
    a3 = 9.6521e-05,
    a4 = -131403,
    c1 = 406.722,
    c2 = 1015448,
    w_ref = 567020,
    z = 2,
    a_0 = 3.72,
    MM = 0.23906) "Eu(But)+2";

  //1293
  constant DataRecord Eu_But2_p1(
    R = Modelica.Constants.R/Eu_But2_p1.MM,
    G_ref = -1309877,
    H_ref = -1705231,
    S_ref = 60.413,
    a1 = 7.4425e-05,
    a2 = 14916.13,
    a3 = -0.00034557,
    a4 = -177933,
    c1 = 824.4882,
    c2 = 2604168,
    w_ref = 140210,
    z = 1,
    a_0 = 3.72,
    MM = 0.32615) "Eu(But)2+";

  //1294
  constant DataRecord Gd_But_p2(
    R = Modelica.Constants.R/Gd_But_p2.MM,
    G_ref = -1032774,
    H_ref = -1236644,
    S_ref = -60.53,
    a1 = 2.8789e-05,
    a2 = 3774.43,
    a3 = 9.2077e-05,
    a4 = -131871,
    c1 = 407.5915,
    c2 = 1028871,
    w_ref = 534550,
    z = 2,
    a_0 = 3.72,
    MM = 0.24435) "Gd(But)+2";

  //1295
  constant DataRecord Gd_But2_p1(
    R = Modelica.Constants.R/Gd_But2_p1.MM,
    G_ref = -1397628,
    H_ref = -1780133,
    S_ref = 91.157,
    a1 = 7.4906e-05,
    a2 = 15032.86,
    a3 = -0.00035002,
    a4 = -178414,
    c1 = 827.6981,
    c2 = 2630355,
    w_ref = 93260,
    z = 1,
    a_0 = 3.72,
    MM = 0.33144) "Gd(But)2+";

  //1296
  constant DataRecord Yb_But_p2(
    R = Modelica.Constants.R/Yb_But_p2.MM,
    G_ref = -1009570,
    H_ref = -1221724,
    S_ref = -106.219,
    a1 = 2.6401e-05,
    a2 = 3191.3,
    a3 = 0.00011498,
    a4 = -129461,
    c1 = 411.0228,
    c2 = 1019281,
    w_ref = 601740,
    z = 2,
    a_0 = 3.72,
    MM = 0.26014) "Yb(But)+2";

  //1297
  constant DataRecord Yb_But2_p1(
    R = Modelica.Constants.R/Yb_But2_p1.MM,
    G_ref = -1374595,
    H_ref = -1773104,
    S_ref = 28.861,
    a1 = 7.2312e-05,
    a2 = 14401.75,
    a3 = -0.00032567,
    a4 = -175807,
    c1 = 831.0524,
    c2 = 2611649,
    w_ref = 188110,
    z = 1,
    a_0 = 3.72,
    MM = 0.34723) "Yb(But)2+";

  //1298
  constant DataRecord Pb_But_p1(
    R = Modelica.Constants.R/Pb_But_p1.MM,
    G_ref = -388346,
    H_ref = -533715,
    S_ref = 197.644,
    a1 = 4.2571e-05,
    a2 = 7139.7,
    a3 = -4.0217e-05,
    a4 = -145783,
    c1 = 368.946,
    c2 = 1087677,
    w_ref = -68580,
    z = 1,
    a_0 = 3.72,
    MM = 0.29429) "Pb(But)+";

  //1299
  constant DataRecord Pb_But2_aq(
    R = Modelica.Constants.R/Pb_But2_aq.MM,
    G_ref = -749258,
    H_ref = -1066414,
    S_ref = 372.301,
    a1 = 9.2175e-05,
    a2 = 19252.09,
    a3 = -0.00051638,
    a4 = -195857,
    c1 = 845.5308,
    c2 = 2726223,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.31738) "Pb(But)2,aq";

  //1300
  constant DataRecord Mn_But_p1(
    R = Modelica.Constants.R/Mn_But_p1.MM,
    G_ref = -589057,
    H_ref = -761647,
    S_ref = 73.283,
    a1 = 4.2252e-05,
    a2 = 7061.34,
    a3 = -3.7016e-05,
    a4 = -145461,
    c1 = 433.7871,
    c2 = 1252593,
    w_ref = 120210,
    z = 1,
    a_0 = 3.72,
    MM = 0.14204) "Mn(But)+";

  //1301
  constant DataRecord Mn_But2_aq(
    R = Modelica.Constants.R/Mn_But2_aq.MM,
    G_ref = -947120,
    H_ref = -1302939,
    S_ref = 209.56,
    a1 = 9.1109e-05,
    a2 = 18990.97,
    a3 = -0.00050595,
    a4 = -194778,
    c1 = 938.1072,
    c2 = 3047990,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16513) "Mn(But)2,aq";

  //1302
  constant DataRecord Ni_But_p1(
    R = Modelica.Constants.R/Ni_But_p1.MM,
    G_ref = -409718,
    H_ref = -604124,
    S_ref = -2.017,
    a1 = 3.5048e-05,
    a2 = 5301.59,
    a3 = 3.2271e-05,
    a4 = -138185,
    c1 = 402.3995,
    c2 = 1106852,
    w_ref = 234640,
    z = 1,
    a_0 = 3.72,
    MM = 0.14581) "Ni(But)+";

  //1303
  constant DataRecord Ni_But2_aq(
    R = Modelica.Constants.R/Ni_But2_aq.MM,
    G_ref = -770346,
    H_ref = -1154910,
    S_ref = 111.018,
    a1 = 8.2651e-05,
    a2 = 16925.16,
    a3 = -0.00042463,
    a4 = -186238,
    c1 = 856.2958,
    c2 = 2763637,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.1689) "Ni(But)2,aq";

  //1304
  constant DataRecord Zn_But_p1(
    R = Modelica.Constants.R/Zn_But_p1.MM,
    G_ref = -509674,
    H_ref = -699761,
    S_ref = 24.225,
    a1 = 3.7906e-05,
    a2 = 5999.27,
    a3 = 4.883e-06,
    a4 = -141072,
    c1 = 428.3621,
    c2 = 1210406,
    w_ref = 193090,
    z = 1,
    a_0 = 3.72,
    MM = 0.15247) "Zn(But)+";

  //1305
  constant DataRecord Zn_But2_aq(
    R = Modelica.Constants.R/Zn_But2_aq.MM,
    G_ref = -868879,
    H_ref = -1246711,
    S_ref = 145.356,
    a1 = 8.5991e-05,
    a2 = 17741.75,
    a3 = -0.00045691,
    a4 = -189615,
    c1 = 914.4249,
    c2 = 2965678,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17556) "Zn(But)2,aq";

  //1306
  constant DataRecord Fe_But_p1(
    R = Modelica.Constants.R/Fe_But_p1.MM,
    G_ref = -454859,
    H_ref = -639646,
    S_ref = 27.648,
    a1 = 3.9228e-05,
    a2 = 6323.74,
    a3 = -8.176e-06,
    a4 = -142411,
    c1 = 417.9724,
    c2 = 1175888,
    w_ref = 188110,
    z = 1,
    a_0 = 3.72,
    MM = 0.14295) "Fe(But)+";

  //1307
  constant DataRecord Fe_But2_aq(
    R = Modelica.Constants.R/Fe_But2_aq.MM,
    G_ref = -814805,
    H_ref = -1187022,
    S_ref = 149.837,
    a1 = 8.7484e-05,
    a2 = 18105.3,
    a3 = -0.000471,
    a4 = -191117,
    c1 = 895.048,
    c2 = 2898332,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16604) "Fe(But)2,aq";

  //1308
  constant DataRecord Co_But_p1(
    R = Modelica.Constants.R/Co_But_p1.MM,
    G_ref = -418333,
    H_ref = -606324,
    S_ref = 19.661,
    a1 = 3.7868e-05,
    a2 = 5990.86,
    a3 = 5.021e-06,
    a4 = -141034,
    c1 = 420.7962,
    c2 = 1181641,
    w_ref = 200790,
    z = 1,
    a_0 = 3.72,
    MM = 0.14603) "Co(But)+";

  //1309
  constant DataRecord Co_But2_aq(
    R = Modelica.Constants.R/Co_But2_aq.MM,
    G_ref = -778789,
    H_ref = -1154943,
    S_ref = 139.386,
    a1 = 8.592e-05,
    a2 = 17723.84,
    a3 = -0.00045609,
    a4 = -189539,
    c1 = 898.278,
    c2 = 2909554,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16912) "Co(But)2,aq";

  //1310
  constant DataRecord Ca_Pent_p1(
    R = Modelica.Constants.R/Ca_Pent_p1.MM,
    G_ref = -900100,
    H_ref = -1101174,
    S_ref = 123.867,
    a1 = 5.0545e-05,
    a2 = 9087.02,
    a3 = -0.00011683,
    a4 = -153833,
    c1 = 563.8337,
    c2 = 1729389,
    w_ref = 42800,
    z = 1,
    a_0 = 3.72,
    MM = 0.1412) "Ca(Pent)+";

  //1311
  constant DataRecord Ca_Pent2_aq(
    R = Modelica.Constants.R/Ca_Pent2_aq.MM,
    G_ref = -1245301,
    H_ref = -1657529,
    S_ref = 302.951,
    a1 = 0.00010981,
    a2 = 23556.63,
    a3 = -0.00068547,
    a4 = -213652,
    c1 = 1273.9435,
    c2 = 4215275,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.17832) "Ca(Pent)2,aq";

  //1312
  constant DataRecord Sr_Pent_p1(
    R = Modelica.Constants.R/Sr_Pent_p1.MM,
    G_ref = -908405,
    H_ref = -1103551,
    S_ref = 157.925,
    a1 = 5.0787e-05,
    a2 = 9143.76,
    a3 = -0.00011856,
    a4 = -154067,
    c1 = 545.2204,
    c2 = 1681064,
    w_ref = -8330,
    z = 1,
    a_0 = 3.72,
    MM = 0.18874) "Sr(Pent)+";

  //1313
  constant DataRecord Sr_Pent2_aq(
    R = Modelica.Constants.R/Sr_Pent2_aq.MM,
    G_ref = -1251321,
    H_ref = -1654487,
    S_ref = 347.519,
    a1 = 0.00011027,
    a2 = 23667.97,
    a3 = -0.00068953,
    a4 = -214112,
    c1 = 1246.8165,
    c2 = 4120989,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.22586) "Sr(Pent)2,aq";

  //1314
  constant DataRecord Ba_Pent_p1(
    R = Modelica.Constants.R/Ba_Pent_p1.MM,
    G_ref = -905238,
    H_ref = -1085715,
    S_ref = 213.999,
    a1 = 5.3568e-05,
    a2 = 9824.95,
    a3 = -0.00014578,
    a4 = -156887,
    c1 = 524.9748,
    c2 = 1637915,
    w_ref = -93350,
    z = 1,
    a_0 = 3.72,
    MM = 0.23846) "Ba(Pent)+";

  //1315
  constant DataRecord Ba_Pent2_aq(
    R = Modelica.Constants.R/Ba_Pent2_aq.MM,
    G_ref = -1248041,
    H_ref = -1631379,
    S_ref = 420.902,
    a1 = 0.00011369,
    a2 = 24503.26,
    a3 = -0.00072247,
    a4 = -217564,
    c1 = 1222.5966,
    c2 = 4036803,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.27558) "Ba(Pent)2,aq";

  //1316
  constant DataRecord Cu_Pent_p1(
    R = Modelica.Constants.R/Cu_Pent_p1.MM,
    G_ref = -292114,
    H_ref = -507189,
    S_ref = 68.534,
    a1 = 4.6658e-05,
    a2 = 8136.71,
    a3 = -7.9216e-05,
    a4 = -149904,
    c1 = 581.7074,
    c2 = 1764481,
    w_ref = 127240,
    z = 1,
    a_0 = 3.72,
    MM = 0.16466) "Cu(Pent)+";

  //1317
  constant DataRecord Cu_Pent2_aq(
    R = Modelica.Constants.R/Cu_Pent2_aq.MM,
    G_ref = -645934,
    H_ref = -1077254,
    S_ref = 230.538,
    a1 = 0.00010516,
    a2 = 22422.14,
    a3 = -0.00064099,
    a4 = -208962,
    c1 = 1293.6434,
    c2 = 4283742,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.20178) "Cu(Pent)2,aq";

  //1318
  constant DataRecord Na_Pent_aq(
    R = Modelica.Constants.R/Na_Pent_aq.MM,
    G_ref = -607764,
    H_ref = -803366,
    S_ref = 217.288,
    a1 = 6.117e-05,
    a2 = 11681.44,
    a3 = -0.00021879,
    a4 = -164561,
    c1 = 539.1841,
    c2 = 1661441,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.092111) "Na(Pent),aq";

  //1319
  constant DataRecord Na_Pent2_n1(
    R = Modelica.Constants.R/Na_Pent2_n1.MM,
    G_ref = -951764,
    H_ref = -1372055,
    S_ref = 350.979,
    a1 = 0.0001224,
    a2 = 26631.95,
    a3 = -0.00080625,
    a4 = -226367,
    c1 = 1245.3433,
    c2 = 4063823,
    w_ref = 149950,
    z = -1,
    a_0 = 3.72,
    MM = 0.22523) "Na(Pent)2-";

  //1320
  constant DataRecord K_Pent_aq(
    R = Modelica.Constants.R/K_Pent_aq.MM,
    G_ref = -628232,
    H_ref = -811198,
    S_ref = 273.081,
    a1 = 6.7658e-05,
    a2 = 13263.49,
    a3 = -0.00028054,
    a4 = -171100,
    c1 = 500.1231,
    c2 = 1525670,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.10822) "K(Pent),aq";

  //1321
  constant DataRecord K_Pent2_n1(
    R = Modelica.Constants.R/K_Pent2_n1.MM,
    G_ref = -972119,
    H_ref = -1375553,
    S_ref = 420.936,
    a1 = 0.00012927,
    a2 = 28310.07,
    a3 = -0.00087239,
    a4 = -233304,
    c1 = 1159.3245,
    c2 = 3798926,
    w_ref = 43560,
    z = -1,
    a_0 = 3.72,
    MM = 0.24134) "K(Pent)2-";

  //1322
  constant DataRecord La_Pent_p2(
    R = Modelica.Constants.R/La_Pent_p2.MM,
    G_ref = -1046787,
    H_ref = -1284793,
    S_ref = -49.949,
    a1 = 3.9046e-05,
    a2 = 6277.55,
    a3 = -6.021e-06,
    a4 = -142223,
    c1 = 554.939,
    c2 = 1546511,
    w_ref = 517390,
    z = 2,
    a_0 = 3.72,
    MM = 0.24003) "La(Pent)+2";

  //1323
  constant DataRecord La_Pent2_p1(
    R = Modelica.Constants.R/La_Pent2_p1.MM,
    G_ref = -1403000,
    H_ref = -1851077,
    S_ref = 132.779,
    a1 = 9.535e-05,
    a2 = 20025.75,
    a3 = -0.00054641,
    a4 = -199054,
    c1 = 1180.6353,
    c2 = 3877334,
    w_ref = 30000,
    z = 1,
    a_0 = 3.72,
    MM = 0.34115) "La(Pent)2+";

  //1324
  constant DataRecord Eu_Pent_p2(
    R = Modelica.Constants.R/Eu_Pent_p2.MM,
    G_ref = -935986,
    H_ref = -1182047,
    S_ref = -55.882,
    a1 = 3.7358e-05,
    a2 = 5866.01,
    a3 = 1.0021e-05,
    a4 = -140520,
    c1 = 559.1916,
    c2 = 1558017,
    w_ref = 527600,
    z = 2,
    a_0 = 3.72,
    MM = 0.25308) "Eu(Pent)+2";

  //1325
  constant DataRecord Eu_Pent2_p1(
    R = Modelica.Constants.R/Eu_Pent2_p1.MM,
    G_ref = -1293002,
    H_ref = -1749774,
    S_ref = 124.692,
    a1 = 9.3471e-05,
    a2 = 19567.69,
    a3 = -0.00052859,
    a4 = -197163,
    c1 = 1188.1794,
    c2 = 3899785,
    w_ref = 41800,
    z = 1,
    a_0 = 3.72,
    MM = 0.3542) "Eu(Pent)2+";

  //1326
  constant DataRecord U_Pent_p2(
    R = Modelica.Constants.R/U_Pent_p2.MM,
    G_ref = -837415,
    H_ref = -1062928,
    S_ref = -14.937,
    a1 = 3.8937e-05,
    a2 = 6251.31,
    a3 = -5.058e-06,
    a4 = -142114,
    c1 = 551.8842,
    c2 = 1552264,
    w_ref = 466260,
    z = 2,
    a_0 = 3.72,
    MM = 0.33915) "U(Pent)+2";

  //1327
  constant DataRecord Eu_Pent_p1(
    R = Modelica.Constants.R/Eu_Pent_p1.MM,
    G_ref = -885811,
    H_ref = -1079003,
    S_ref = 187.188,
    a1 = 6.034e-05,
    a2 = 11478.8,
    a3 = -0.00021084,
    a4 = -163724,
    c1 = 629.692,
    c2 = 1988843,
    w_ref = -52590,
    z = 1,
    a_0 = 3.72,
    MM = 0.25308) "Eu(Pent)+";

  //1328
  constant DataRecord Gd_Pent_p2(
    R = Modelica.Constants.R/Gd_Pent_p2.MM,
    G_ref = -1024365,
    H_ref = -1260907,
    S_ref = -33.334,
    a1 = 3.782e-05,
    a2 = 5980.07,
    a3 = 5.268e-06,
    a4 = -140992,
    c1 = 559.9878,
    c2 = 1571439,
    w_ref = 494340,
    z = 2,
    a_0 = 3.72,
    MM = 0.25837) "Gd(Pent)+2";

  //1329
  constant DataRecord Gd_Pent2_p1(
    R = Modelica.Constants.R/Gd_Pent2_p1.MM,
    G_ref = -1380749,
    H_ref = -1825563,
    S_ref = 155.436,
    a1 = 9.3953e-05,
    a2 = 19684.21,
    a3 = -0.00053291,
    a4 = -197644,
    c1 = 1191.4174,
    c2 = 3925973,
    w_ref = -4850,
    z = 1,
    a_0 = 3.72,
    MM = 0.35949) "Gd(Pent)2+";

  //1330
  constant DataRecord Yb_Pent_p2(
    R = Modelica.Constants.R/Yb_Pent_p2.MM,
    G_ref = -1001160,
    H_ref = -1248836,
    S_ref = -79.023,
    a1 = 3.5437e-05,
    a2 = 5398.03,
    a3 = 2.8221e-05,
    a4 = -138587,
    c1 = 563.5852,
    c2 = 1561850,
    w_ref = 563330,
    z = 2,
    a_0 = 3.72,
    MM = 0.27416) "Yb(Pent)+2";

  //1331
  constant DataRecord Yb_Pent2_p1(
    R = Modelica.Constants.R/Yb_Pent2_p1.MM,
    G_ref = -1357721,
    H_ref = -1818613,
    S_ref = 93.136,
    a1 = 9.1361e-05,
    a2 = 19052.81,
    a3 = -0.00050842,
    a4 = -195033,
    c1 = 1194.8081,
    c2 = 3907266,
    w_ref = 90370,
    z = 1,
    a_0 = 3.72,
    MM = 0.37528) "Yb(Pent)2+";

  //1332
  constant DataRecord Mg_Pent_p1(
    R = Modelica.Constants.R/Mg_Pent_p1.MM,
    G_ref = -801353,
    H_ref = -1032946,
    S_ref = 12.627,
    a1 = 4.8888e-05,
    a2 = 8680.34,
    a3 = -0.00010033,
    a4 = -152155,
    c1 = 591.4574,
    c2 = 1771384,
    w_ref = 211500,
    z = 1,
    a_0 = 3.72,
    MM = 0.12543) "Mg(Pent)+";

  //1333
  constant DataRecord Mg_Pent2_aq(
    R = Modelica.Constants.R/Mg_Pent2_aq.MM,
    G_ref = -1146608,
    H_ref = -1599598,
    S_ref = 157.381,
    a1 = 0.00010732,
    a2 = 22950.75,
    a3 = -0.00066161,
    a4 = -211146,
    c1 = 1297.5182,
    c2 = 4297215,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.16255) "Mg(Pent)2,aq";

  //1334
  constant DataRecord Pb_Pent_p1(
    R = Modelica.Constants.R/Pb_Pent_p1.MM,
    G_ref = -379878,
    H_ref = -557819,
    S_ref = 224.84,
    a1 = 5.1599e-05,
    a2 = 9343.5,
    a3 = -0.0001267,
    a4 = -154896,
    c1 = 521.2674,
    c2 = 1630245,
    w_ref = -109620,
    z = 1,
    a_0 = 3.72,
    MM = 0.30831) "Pb(Pent)+";

  //1335
  constant DataRecord Pb_Pent2_aq(
    R = Modelica.Constants.R/Pb_Pent2_aq.MM,
    G_ref = -732267,
    H_ref = -1112078,
    S_ref = 435.086,
    a1 = 0.00011155,
    a2 = 23983.61,
    a3 = -0.00070226,
    a4 = -215417,
    c1 = 1218.2908,
    c2 = 4021837,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.34543) "Pb(Pent)2,aq";

  //1336
  constant DataRecord Co_Pent_p1(
    R = Modelica.Constants.R/Co_Pent_p1.MM,
    G_ref = -409752,
    H_ref = -630416,
    S_ref = 46.857,
    a1 = 4.6898e-05,
    a2 = 8196.5,
    a3 = -8.1801e-05,
    a4 = -150155,
    c1 = 573.1896,
    c2 = 1724210,
    w_ref = 160540,
    z = 1,
    a_0 = 3.72,
    MM = 0.16005) "Co(Pent)+";

  //1337
  constant DataRecord Co_Pent2_aq(
    R = Modelica.Constants.R/Co_Pent2_aq.MM,
    G_ref = -761630,
    H_ref = -1200536,
    S_ref = 202.171,
    a1 = 0.0001053,
    a2 = 22455.36,
    a3 = -0.00064198,
    a4 = -209100,
    c1 = 1271.0377,
    c2 = 4205171,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19717) "Co(Pent)2,aq";

  //1338
  constant DataRecord Ni_Pent_p1(
    R = Modelica.Constants.R/Ni_Pent_p1.MM,
    G_ref = -401137,
    H_ref = -628127,
    S_ref = 25.179,
    a1 = 4.4074e-05,
    a2 = 7505.8,
    a3 = -5.4421e-05,
    a4 = -147298,
    c1 = 554.6712,
    c2 = 1649421,
    w_ref = 193090,
    z = 1,
    a_0 = 3.72,
    MM = 0.15983) "Ni(Pent)+";

  //1339
  constant DataRecord Ni_Pent2_aq(
    R = Modelica.Constants.R/Ni_Pent2_aq.MM,
    G_ref = -753128,
    H_ref = -1200356,
    S_ref = 173.803,
    a1 = 0.00010203,
    a2 = 21656.68,
    a3 = -0.00061052,
    a4 = -205798,
    c1 = 1229.055,
    c2 = 4059254,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19695) "Ni(Pent)2,aq";

  //1340
  constant DataRecord Zn_Pent_p1(
    R = Modelica.Constants.R/Zn_Pent_p1.MM,
    G_ref = -500754,
    H_ref = -723397,
    S_ref = 51.421,
    a1 = 4.6934e-05,
    a2 = 8202.94,
    a3 = -8.1563e-05,
    a4 = -150180,
    c1 = 580.6923,
    c2 = 1752975,
    w_ref = 152130,
    z = 1,
    a_0 = 3.72,
    MM = 0.16649) "Zn(Pent)+";

  //1341
  constant DataRecord Zn_Pent2_aq(
    R = Modelica.Constants.R/Zn_Pent2_aq.MM,
    G_ref = -851088,
    H_ref = -1291559,
    S_ref = 208.141,
    a1 = 0.00010537,
    a2 = 22473.31,
    a3 = -0.00064281,
    a4 = -209175,
    c1 = 1287.1842,
    c2 = 4261295,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.20361) "Zn(Pent)2,aq";

  //1342
  constant DataRecord Cd_Pent_p1(
    R = Modelica.Constants.R/Cd_Pent_p1.MM,
    G_ref = -432559,
    H_ref = -643349,
    S_ref = 101.621,
    a1 = 5.2228e-05,
    a2 = 9497.05,
    a3 = -0.00013271,
    a4 = -155532,
    c1 = 589.1842,
    c2 = 1806668,
    w_ref = 76650,
    z = 1,
    a_0 = 3.72,
    MM = 0.21352) "Cd(Pent)+";

  //1343
  constant DataRecord Cd_Pent2_aq(
    R = Modelica.Constants.R/Cd_Pent2_aq.MM,
    G_ref = -784036,
    H_ref = -1208030,
    S_ref = 273.839,
    a1 = 0.00011155,
    a2 = 23983.61,
    a3 = -0.00070226,
    a4 = -215417,
    c1 = 1317.3249,
    c2 = 4366058,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.25064) "Cd(Pent)2,aq";

  //1344
  constant DataRecord Fe_Pent_p1(
    R = Modelica.Constants.R/Fe_Pent_p1.MM,
    G_ref = -446161,
    H_ref = -663515,
    S_ref = 54.844,
    a1 = 4.826e-05,
    a2 = 8529.17,
    a3 = -9.4922e-05,
    a4 = -151528,
    c1 = 570.3876,
    c2 = 1718457,
    w_ref = 148070,
    z = 1,
    a_0 = 3.72,
    MM = 0.15697) "Fe(Pent)+";

  //1345
  constant DataRecord Fe_Pent2_aq(
    R = Modelica.Constants.R/Fe_Pent2_aq.MM,
    G_ref = -797412,
    H_ref = -1230054,
    S_ref = 212.547,
    a1 = 0.00010686,
    a2 = 22836.82,
    a3 = -0.00065689,
    a4 = -210677,
    c1 = 1267.8076,
    c2 = 4193950,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19409) "Fe(Pent)2,aq";

  //1346
  constant DataRecord Mn_Pent_p1(
    R = Modelica.Constants.R/Mn_Pent_p1.MM,
    G_ref = -579907,
    H_ref = -785111,
    S_ref = 100.479,
    a1 = 5.128e-05,
    a2 = 9265.05,
    a3 = -0.00012347,
    a4 = -154570,
    c1 = 586.1177,
    c2 = 1795162,
    w_ref = 79290,
    z = 1,
    a_0 = 3.72,
    MM = 0.15606) "Mn(Pent)+";

  //1347
  constant DataRecord Mn_Pent2_aq(
    R = Modelica.Constants.R/Mn_Pent2_aq.MM,
    G_ref = -928873,
    H_ref = -1347386,
    S_ref = 272.345,
    a1 = 0.00011049,
    a2 = 23722.49,
    a3 = -0.00069183,
    a4 = -214338,
    c1 = 1310.8669,
    c2 = 4343607,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.19318) "Mn(Pent)2,aq";

  //1348
  constant DataRecord NaSO4_n1(
    R = Modelica.Constants.R/NaSO4_n1.MM,
    G_ref = -1010336,
    H_ref = -1145299,
    S_ref = 106.525,
    a1 = 2.7044e-05,
    a2 = 3349.08,
    a3 = 0.00010866,
    a4 = -130114,
    c1 = 60.6379,
    c2 = -172305,
    w_ref = 519740,
    z = -1,
    a_0 = 3.72,
    MM = 0.11905) "NaSO4-";

  //1349
  constant DataRecord MgSO4_aq(
    R = Modelica.Constants.R/MgSO4_aq.MM,
    G_ref = -1211172,
    H_ref = -1369704,
    S_ref = -56.484,
    a1 = 8.104e-06,
    a2 = -1278.09,
    a3 = 0.00029101,
    a4 = -110985,
    c1 = -28.6144,
    c2 = -312080,
    w_ref = -12550,
    z = 0,
    a_0 = 3.72,
    MM = 0.12037) "MgSO4,aq";

  //1350
  constant DataRecord HCl_aq(
    R = Modelica.Constants.R/HCl_aq.MM,
    G_ref = -127235,
    H_ref = -175954,
    S_ref = 13.389,
    a1 = 5.2497e-06,
    a2 = -1973.89,
    a3 = 0.00031816,
    a4 = -108115,
    c1 = 69.9289,
    c2 = 120194,
    w_ref = -292880,
    z = 0,
    a_0 = 3.72,
    MM = 0.036458) "HCl,aq";

  //1351
  constant DataRecord MgCO3_aq(
    R = Modelica.Constants.R/MgCO3_aq.MM,
    G_ref = -998972,
    H_ref = -1132065,
    S_ref = -100.416,
    a1 = -3.0773e-06,
    a2 = -4005.97,
    a3 = 0.00039774,
    a4 = -99709,
    c1 = -42.8509,
    c2 = -360489,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.08432) "MgCO3,aq";

  //1352
  constant DataRecord CaCO3_aq(
    R = Modelica.Constants.R/CaCO3_aq.MM,
    G_ref = -1099764,
    H_ref = -1202440,
    S_ref = 10.46,
    a1 = -1.6347e-06,
    a2 = -3653.68,
    a3 = 0.0003839,
    a4 = -101165,
    c1 = -48.2453,
    c2 = -379242,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.10009) "CaCO3,aq";

  //1353
  constant DataRecord SrCO3_aq(
    R = Modelica.Constants.R/SrCO3_aq.MM,
    G_ref = -1108174,
    H_ref = -1207586,
    S_ref = 35.564,
    a1 = -1.3941e-06,
    a2 = -3594.98,
    a3 = 0.00038158,
    a4 = -101408,
    c1 = -54.3757,
    c2 = -400547,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.14763) "SrCO3,aq";

  //1354
  constant DataRecord BaCO3_aq(
    R = Modelica.Constants.R/BaCO3_aq.MM,
    G_ref = -1103865,
    H_ref = -1195996,
    S_ref = 66.944,
    a1 = 5.694e-07,
    a2 = -3115.45,
    a3 = 0.00036274,
    a4 = -103391,
    c1 = -60.0153,
    c2 = -420149,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.19735) "BaCO3,aq";

  //1355
  constant DataRecord BeCl_p1(
    R = Modelica.Constants.R/BeCl_p1.MM,
    G_ref = -482482,
    H_ref = -522992,
    S_ref = 386.225,
    a1 = 5.0099e-06,
    a2 = -2031.25,
    a3 = 0.00032013,
    a4 = -107876,
    c1 = 69.6531,
    c2 = 138896,
    w_ref = -354260,
    z = 1,
    a_0 = 3.72,
    MM = 0.04446) "BeCl+";

  //1356
  constant DataRecord BeCl2_aq(
    R = Modelica.Constants.R/BeCl2_aq.MM,
    G_ref = -608860,
    H_ref = -654511,
    S_ref = 545.761,
    a1 = 2.1044e-05,
    a2 = 1883.89,
    a3 = 0.00016624,
    a4 = -124060,
    c1 = 165.9082,
    c2 = 365104,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.07991) "BeCl2,aq";

  //1357
  constant DataRecord FeCl_p2(
    R = Modelica.Constants.R/FeCl_p2.MM,
    G_ref = -156975,
    H_ref = -212631,
    S_ref = -178.824,
    a1 = -2.9974e-06,
    a2 = -3986.39,
    a3 = 0.00039697,
    a4 = -99793,
    c1 = 99.6415,
    c2 = -98249,
    w_ref = 711820,
    z = 2,
    a_0 = 3.72,
    MM = 0.0913) "FeCl+2";

  //1358
  constant DataRecord CoCl_p1(
    R = Modelica.Constants.R/CoCl_p1.MM,
    G_ref = -188937,
    H_ref = -225790,
    S_ref = -47.154,
    a1 = 7.5429e-06,
    a2 = -1412.77,
    a3 = 0.00029582,
    a4 = -110432,
    c1 = 95.2513,
    c2 = 18087,
    w_ref = 300870,
    z = 1,
    a_0 = 3.72,
    MM = 0.09438) "CoCl+";

  //1359
  constant DataRecord CuCl_aq(
    R = Modelica.Constants.R/CuCl_aq.MM,
    G_ref = -94592,
    H_ref = -110198,
    S_ref = 92.299,
    a1 = 1.7189e-05,
    a2 = 942.66,
    a3 = 0.00020324,
    a4 = -120169,
    c1 = 72.5054,
    c2 = 40459,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.09899) "CuCl,aq";

  //1360
  constant DataRecord CuCl2_n1(
    R = Modelica.Constants.R/CuCl2_n1.MM,
    G_ref = -242831,
    H_ref = -305026,
    S_ref = 112.801,
    a1 = 3.5122e-05,
    a2 = 5321.29,
    a3 = 3.1137e-05,
    a4 = -138269,
    c1 = 153.785,
    c2 = 154164,
    w_ref = 511240,
    z = -1,
    a_0 = 3.72,
    MM = 0.13444) "CuCl2-";

  //1361
  constant DataRecord CuCl3_n2(
    R = Modelica.Constants.R/CuCl3_n2.MM,
    G_ref = -376405,
    H_ref = -495904,
    S_ref = 97.404,
    a1 = 5.543e-05,
    a2 = 10280.09,
    a3 = -0.00016376,
    a4 = -158770,
    c1 = 267.6957,
    c2 = 330908,
    w_ref = 1195750,
    z = -2,
    a_0 = 3.72,
    MM = 0.16989) "CuCl3-2";

  //1362
  constant DataRecord CuCl_p1(
    R = Modelica.Constants.R/CuCl_p1.MM,
    G_ref = -67990,
    H_ref = -99776,
    S_ref = -27.238,
    a1 = 7.3136e-06,
    a2 = -1468.71,
    a3 = 0.00029802,
    a4 = -110198,
    c1 = 122.8255,
    c2 = 123554,
    w_ref = 270790,
    z = 1,
    a_0 = 3.72,
    MM = 0.09899) "CuCl+";

  //1363
  constant DataRecord CuCl2_aq(
    R = Modelica.Constants.R/CuCl2_aq.MM,
    G_ref = -193058,
    H_ref = -268445,
    S_ref = 3.264,
    a1 = 2.1257e-05,
    a2 = 1935.94,
    a3 = 0.0001642,
    a4 = -124273,
    c1 = 157.2966,
    c2 = 335172,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.13444) "CuCl2,aq";

  //1364
  constant DataRecord CuCl3_n1(
    R = Modelica.Constants.R/CuCl3_n1.MM,
    G_ref = -315214,
    H_ref = -447044,
    S_ref = -9.247,
    a1 = 4.0275e-05,
    a2 = 6579.42,
    a3 = -1.8313e-05,
    a4 = -143469,
    c1 = 267.3078,
    c2 = 489980,
    w_ref = 694750,
    z = -1,
    a_0 = 3.72,
    MM = 0.16989) "CuCl3-";

  //1365
  constant DataRecord CuCl4_n2(
    R = Modelica.Constants.R/CuCl4_n2.MM,
    G_ref = -433375,
    H_ref = -638202,
    S_ref = -77.32,
    a1 = 6.1378e-05,
    a2 = 11732.23,
    a3 = -0.00022084,
    a4 = -164774,
    c1 = 366.1155,
    c2 = 587986,
    w_ref = 1461180,
    z = -2,
    a_0 = 3.72,
    MM = 0.20534) "CuCl4-2";

  //1366
  constant DataRecord CdCl_p1(
    R = Modelica.Constants.R/CdCl_p1.MM,
    G_ref = -220191,
    H_ref = -248496,
    S_ref = 3.138,
    a1 = 1.2966e-05,
    a2 = -88.53,
    a3 = 0.00024377,
    a4 = -115905,
    c1 = 111.4973,
    c2 = 98625,
    w_ref = 225680,
    z = 1,
    a_0 = 4.25,
    MM = 0.14785) "CdCl+";

  //1367
  constant DataRecord CdCl2_aq(
    R = Modelica.Constants.R/CdCl2_aq.MM,
    G_ref = -355017,
    H_ref = -424099,
    S_ref = 43.137,
    a1 = 2.7725e-05,
    a2 = 3515.23,
    a3 = 0.00010212,
    a4 = -130804,
    c1 = 143.3024,
    c2 = 286529,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.1833) "CdCl2,aq";

  //1368
  constant DataRecord CdCl3_n1(
    R = Modelica.Constants.R/CdCl3_n1.MM,
    G_ref = -485223,
    H_ref = -606329,
    S_ref = 45.438,
    a1 = 4.7206e-05,
    a2 = 8271.98,
    a3 = -8.4835e-05,
    a4 = -150469,
    c1 = 236.1663,
    c2 = 407768,
    w_ref = 613460,
    z = -1,
    a_0 = 3.72,
    MM = 0.21875) "CdCl3-";

  //1369
  constant DataRecord CdCl4_n2(
    R = Modelica.Constants.R/CdCl4_n2.MM,
    G_ref = -611203,
    H_ref = -798278,
    S_ref = 0.962,
    a1 = 6.9002e-05,
    a2 = 13593.94,
    a3 = -0.00029401,
    a4 = -172469,
    c1 = 318.9279,
    c2 = 462340,
    w_ref = 1341390,
    z = -2,
    a_0 = 3.72,
    MM = 0.2542) "CdCl4-2";

  //1370
  constant DataRecord TlCl_aq(
    R = Modelica.Constants.R/TlCl_aq.MM,
    G_ref = -166586,
    H_ref = -158235,
    S_ref = 203.719,
    a1 = 3.3967e-05,
    a2 = 5039.21,
    a3 = 4.2225e-05,
    a4 = -137101,
    c1 = -26.2521,
    c2 = -302796,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.23982) "TlCl,aq";

  //1371
  constant DataRecord TlCl_p2(
    R = Modelica.Constants.R/TlCl_p2.MM,
    G_ref = 39225,
    H_ref = 2807,
    S_ref = -77.488,
    a1 = -1.215e-06,
    a2 = -3551.25,
    a3 = 0.00037987,
    a4 = -101592,
    c1 = 56.3786,
    c2 = -199886,
    w_ref = 559610,
    z = 0,
    a_0 = 3.72,
    MM = 0.23982) "TlCl+2";

  //1372
  constant DataRecord AuCl_aq(
    R = Modelica.Constants.R/AuCl_aq.MM,
    G_ref = -13322,
    H_ref = -8954,
    S_ref = 173.51,
    a1 = 2.9437e-05,
    a2 = 3933.29,
    a3 = 8.5693e-05,
    a4 = -132532,
    c1 = 0.2305,
    c2 = -210752,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.23242) "AuCl,aq";

  //1373
  constant DataRecord AuCl2_n1(
    R = Modelica.Constants.R/AuCl2_n1.MM,
    G_ref = -153892,
    H_ref = -187129,
    S_ref = 224.262,
    a1 = 4.8196e-05,
    a2 = 8513.69,
    a3 = -9.4337e-05,
    a4 = -151465,
    c1 = -2.8301,
    c2 = -335975,
    w_ref = 341960,
    z = -1,
    a_0 = 3.72,
    MM = 0.26787) "AuCl2-";

  //1374
  constant DataRecord AuCl3_n2(
    R = Modelica.Constants.R/AuCl3_n2.MM,
    G_ref = -283721,
    H_ref = -359920,
    S_ref = 256.856,
    a1 = 6.9822e-05,
    a2 = 13794.1,
    a3 = -0.00030188,
    a4 = -173297,
    c1 = 7.1777,
    c2 = -497519,
    w_ref = 955080,
    z = -2,
    a_0 = 3.72,
    MM = 0.30332) "AuCl3-2";

  //1375
  constant DataRecord AuCl4_n1(
    R = Modelica.Constants.R/AuCl4_n1.MM,
    G_ref = -147829,
    H_ref = 98629,
    S_ref = -136.9,
    a1 = -2.243e-06,
    a2 = -3802.25,
    a3 = 0.00038973,
    a4 = -100554,
    c1 = 82.5298,
    c2 = -138520,
    w_ref = 651830,
    z = -1,
    a_0 = 3.72,
    MM = 0.33877) "AuCl4-";

  //1376
  constant DataRecord HgCl_p1(
    R = Modelica.Constants.R/HgCl_p1.MM,
    G_ref = -9017,
    H_ref = -30853,
    S_ref = 48.827,
    a1 = 1.0435e-05,
    a2 = -706.51,
    a3 = 0.00026806,
    a4 = -113349,
    c1 = 136.0001,
    c2 = 206012,
    w_ref = 156270,
    z = 1,
    a_0 = 3.72,
    MM = 0.23604) "HgCl+";

  //1377
  constant DataRecord HgCl2_aq(
    R = Modelica.Constants.R/HgCl2_aq.MM,
    G_ref = -179234,
    H_ref = -237597,
    S_ref = 103.094,
    a1 = 2.5167e-05,
    a2 = 2890.47,
    a3 = 0.00012668,
    a4 = -128219,
    c1 = 203.5424,
    c2 = 496055,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.27149) "HgCl2,aq";

  //1378
  constant DataRecord HgCl3_n1(
    R = Modelica.Constants.R/HgCl3_n1.MM,
    G_ref = -311725,
    H_ref = -415467,
    S_ref = 127.654,
    a1 = 4.3933e-05,
    a2 = 7472.75,
    a3 = -5.3425e-05,
    a4 = -147164,
    c1 = 326.5436,
    c2 = 761911,
    w_ref = 488520,
    z = -1,
    a_0 = 3.72,
    MM = 0.30694) "HgCl3-";

  //1379
  constant DataRecord HgCl4_n2(
    R = Modelica.Constants.R/HgCl4_n2.MM,
    G_ref = -447638,
    H_ref = -606776,
    S_ref = 118.658,
    a1 = 6.523e-05,
    a2 = 12672.83,
    a3 = -0.00025781,
    a4 = -168661,
    c1 = 458.4032,
    c2 = 1003574,
    w_ref = 1165080,
    z = -2,
    a_0 = 3.72,
    MM = 0.34239) "HgCl4-2";

  //1380
  constant DataRecord InCl_p2(
    R = Modelica.Constants.R/InCl_p2.MM,
    G_ref = -247802,
    H_ref = -277018,
    S_ref = -162.339,
    a1 = -2.6949e-06,
    a2 = -3912.58,
    a3 = 0.00039407,
    a4 = -100098,
    c1 = 93.0329,
    c2 = -113591,
    w_ref = 688020,
    z = 2,
    a_0 = 3.72,
    MM = 0.15027) "InCl+2";

  //1381
  constant DataRecord BeF_p1(
    R = Modelica.Constants.R/BeF_p1.MM,
    G_ref = -663135,
    H_ref = -741915,
    S_ref = 247.693,
    a1 = -5.2731e-06,
    a2 = -4542.15,
    a3 = 0.00041881,
    a4 = -97496,
    c1 = 99.4704,
    c2 = 175280,
    w_ref = -144220,
    z = 1,
    a_0 = 3.72,
    MM = 0.02801) "BeF+";

  //1382
  constant DataRecord BeF2_aq(
    R = Modelica.Constants.R/BeF2_aq.MM,
    G_ref = -968119,
    H_ref = -1080936,
    S_ref = 300.076,
    a1 = -2.1966e-06,
    a2 = -3790.87,
    a3 = 0.00038929,
    a4 = -100600,
    c1 = 190.9038,
    c2 = 451981,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.04701) "BeF2,aq";

  //1383
  constant DataRecord BeF3_1(
    R = Modelica.Constants.R/BeF3_1.MM,
    G_ref = -1265794,
    H_ref = -1387046,
    S_ref = 438.4,
    a1 = 8.586e-07,
    a2 = -3044.95,
    a3 = 0.00035997,
    a4 = -103684,
    c1 = 262.9874,
    c2 = 692059,
    w_ref = 16780,
    z = -1,
    a_0 = 3.72,
    MM = 0.06601) "BeF3-1";

  //1384
  constant DataRecord BeF4_n2(
    R = Modelica.Constants.R/BeF4_n2.MM,
    G_ref = -1550971,
    H_ref = -1615426,
    S_ref = 795.462,
    a1 = 4.553e-06,
    a2 = -2142.79,
    a3 = 0.00032451,
    a4 = -107412,
    c1 = 332.8087,
    c2 = 895514,
    w_ref = 139290,
    z = -2,
    a_0 = 3.72,
    MM = 0.08501) "BeF4-2";

  //1385
  constant DataRecord MnF_p1(
    R = Modelica.Constants.R/MnF_p1.MM,
    G_ref = -517314,
    H_ref = -554188,
    S_ref = -55.647,
    a1 = 1.3121e-06,
    a2 = -2934.16,
    a3 = 0.00035562,
    a4 = -104140,
    c1 = 126.7083,
    c2 = 123503,
    w_ref = 313090,
    z = 1,
    a_0 = 3.72,
    MM = 0.07394) "MnF+";

  //1386
  constant DataRecord FeF_p1(
    R = Modelica.Constants.R/FeF_p1.MM,
    G_ref = -381401,
    H_ref = -424475,
    S_ref = -87.529,
    a1 = -1.7715e-06,
    a2 = -3687.11,
    a3 = 0.00038521,
    a4 = -101027,
    c1 = 110.3676,
    c2 = 50635,
    w_ref = 363300,
    z = 1,
    a_0 = 3.72,
    MM = 0.07485) "FeF+";

  //1387
  constant DataRecord FeF_p2(
    R = Modelica.Constants.R/FeF_p2.MM,
    G_ref = -333235,
    H_ref = -364598,
    S_ref = -107.529,
    a1 = -1.4349e-05,
    a2 = -6758,
    a3 = 0.00050591,
    a4 = -88333,
    c1 = 100.3344,
    c2 = -61869,
    w_ref = 605720,
    z = 2,
    a_0 = 3.72,
    MM = 0.07485) "FeF+2";

  //1388
  constant DataRecord CoF_p1(
    R = Modelica.Constants.R/CoF_p1.MM,
    G_ref = -341967,
    H_ref = -389911,
    S_ref = -94.558,
    a1 = -3.2045e-06,
    a2 = -4036.97,
    a3 = 0.00039896,
    a4 = -99583,
    c1 = 112.4061,
    c2 = -54467,
    w_ref = 373460,
    z = 1,
    a_0 = 3.72,
    MM = 0.07793) "CoF+";

  //1389
  constant DataRecord NiF_p1(
    R = Modelica.Constants.R/NiF_p1.MM,
    G_ref = -333749,
    H_ref = -386246,
    S_ref = -110.29,
    a1 = -6.048e-06,
    a2 = -4731.27,
    a3 = 0.00042625,
    a4 = -96713,
    c1 = 93.3743,
    c2 = -20322,
    w_ref = 400410,
    z = 1,
    a_0 = 3.72,
    MM = 0.07771) "NiF+";

  //1390
  constant DataRecord CuF_p1(
    R = Modelica.Constants.R/CuF_p1.MM,
    G_ref = -224844,
    H_ref = -268968,
    S_ref = -78.827,
    a1 = -3.4158e-06,
    a2 = -4088.56,
    a3 = 0.00040099,
    a4 = -99370,
    c1 = 140.4707,
    c2 = 159938,
    w_ref = 348690,
    z = 1,
    a_0 = 3.72,
    MM = 0.08254) "CuF+";

  //1391
  constant DataRecord ZnF_p1(
    R = Modelica.Constants.R/ZnF_p1.MM,
    G_ref = -435592,
    H_ref = -485892,
    S_ref = -91.253,
    a1 = -3.0941e-06,
    a2 = -4010.07,
    a3 = 0.0003979,
    a4 = -99692,
    c1 = 119.6586,
    c2 = 81316,
    w_ref = 368320,
    z = 1,
    a_0 = 3.72,
    MM = 0.08437) "ZnF+";

  //1392
  constant DataRecord AgF_aq(
    R = Modelica.Constants.R/AgF_aq.MM,
    G_ref = -206936,
    H_ref = -228957,
    S_ref = 70.04,
    a1 = 1.079e-05,
    a2 = -619.99,
    a3 = 0.00026466,
    a4 = -113709,
    c1 = 50.9733,
    c2 = -34384,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.12687) "AgF,aq";

  //1393
  constant DataRecord CdF_p1(
    R = Modelica.Constants.R/CdF_p1.MM,
    G_ref = -365569,
    H_ref = -408124,
    S_ref = -54.852,
    a1 = 2.269e-06,
    a2 = -2700.52,
    a3 = 0.00034643,
    a4 = -105106,
    c1 = 130.0186,
    c2 = 135009,
    w_ref = 313090,
    z = 1,
    a_0 = 3.72,
    MM = 0.1314) "CdF+";

  //1394
  constant DataRecord CdF2_aq(
    R = Modelica.Constants.R/CdF2_aq.MM,
    G_ref = -649206,
    H_ref = -754664,
    S_ref = -99.244,
    a1 = 4.4844e-06,
    a2 = -2159.53,
    a3 = 0.00032517,
    a4 = -107345,
    c1 = 168.2981,
    c2 = 373409,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.1504) "CdF2,aq";

  //1395
  // constant DataRecord BaF_p1(
  //   R = Modelica.Constants.R/BaF_p1.MM,
  //   G_ref = -844360,
  //   H_ref = -865799,
  //   S_ref = 26.694,
  //   a1 = 3.8325e-06,
  //   a2 = -2318.77,
  //   a3 = 0.00033143,
  //   a4 = -106684,
  //   c1 = 70.1787,
  //   c2 = -33744,
  //   w_ref = 190580,
  //   z = 1,
  //   a_0 = 3.72,
  //   MM = 0.15634) "BaF+";

  //1396
  constant DataRecord TlF_aq(
    R = Modelica.Constants.R/TlF_aq.MM,
    G_ref = -314704,
    H_ref = -322662,
    S_ref = 138.825,
    a1 = 2.2974e-05,
    a2 = 2355.05,
    a3 = 0.00014772,
    a4 = -126005,
    c1 = -15.7846,
    c2 = -266416,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.22337) "TlF,aq";

  //1397
  constant DataRecord HgF_p1(
    R = Modelica.Constants.R/HgF_p1.MM,
    G_ref = -126374,
    H_ref = -165318,
    S_ref = -18.744,
    a1 = -2.059e-07,
    a2 = -3304.77,
    a3 = 0.00037018,
    a4 = -102608,
    c1 = 156.0578,
    c2 = 242396,
    w_ref = 260370,
    z = 1,
    a_0 = 3.72,
    MM = 0.21959) "HgF+";

  //1398
  constant DataRecord InF_p2(
    R = Modelica.Constants.R/InF_p2.MM,
    G_ref = -405915,
    H_ref = -413036,
    S_ref = -98.45,
    a1 = -1.4019e-05,
    a2 = -6677.66,
    a3 = 0.00050275,
    a4 = -88667,
    c1 = 94.4626,
    c2 = -77207,
    w_ref = 589900,
    z = 2,
    a_0 = 3.72,
    MM = 0.13382) "InF+2";

  //1399
  constant DataRecord PbF_p1(
    R = Modelica.Constants.R/PbF_p1.MM,
    G_ref = -317398,
    H_ref = -337193,
    S_ref = 34.56,
    a1 = 1.8778e-06,
    a2 = -2796.04,
    a3 = 0.00035018,
    a4 = -104713,
    c1 = 66.8578,
    c2 = -41413,
    w_ref = 178490,
    z = 1,
    a_0 = 3.72,
    MM = 0.22619) "PbF+";

  //1400
  // constant DataRecord PbF_p1(
  //   R = Modelica.Constants.R/PbF_p1.MM,
  //   G_ref = -317398,
  //   H_ref = -337193,
  //   S_ref = 34.56,
  //   a1 = 1.8778e-06,
  //   a2 = -2796.04,
  //   a3 = 0.00035018,
  //   a4 = -104713,
  //   c1 = 66.8578,
  //   c2 = -41413,
  //   w_ref = 178490,
  //   z = 1,
  //   a_0 = 3.72,
  //   MM = 0.22619) "PbF+";

  //1401
  constant DataRecord PbF2_aq(
    R = Modelica.Constants.R/PbF2_aq.MM,
    G_ref = -606914,
    H_ref = -681084,
    S_ref = 18.744,
    a1 = 4.5555e-06,
    a2 = -2142.17,
    a3 = 0.00032449,
    a4 = -107416,
    c1 = 69.2632,
    c2 = 29188,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.24519) "PbF2,aq";

  //1402
  constant DataRecord Ag_HS2_n1(
    R = Modelica.Constants.R/Ag_HS2_n1.MM,
    G_ref = 0,
    H_ref = -34309,
    S_ref = 186.941,
    a1 = 4.3369e-05,
    a2 = 7335.01,
    a3 = -4.8011e-05,
    a4 = -146595,
    c1 = 155.9599,
    c2 = 197790,
    w_ref = 398610,
    z = -1,
    a_0 = 3.72,
    MM = 0.17401) "Ag(HS)2-";

  //1403
  constant DataRecord Au_HS2_n1(
    R = Modelica.Constants.R/Au_HS2_n1.MM,
    G_ref = 10163,
    H_ref = -10498,
    S_ref = 237.526,
    a1 = 5.1639e-05,
    a2 = 9354.25,
    a3 = -0.00012737,
    a4 = -154942,
    c1 = 70.3071,
    c2 = -75341,
    w_ref = 321880,
    z = -1,
    a_0 = 3.72,
    MM = 0.26311) "Au(HS)2-";

  //1404
  constant DataRecord Pb_HS2_aq(
    R = Modelica.Constants.R/Pb_HS2_aq.MM,
    G_ref = -83923,
    H_ref = -95609,
    S_ref = 219.911,
    a1 = 3.1307e-05,
    a2 = 4389.73,
    a3 = 6.7752e-05,
    a4 = -134419,
    c1 = 119.2545,
    c2 = 202945,
    w_ref = -15900,
    z = 0,
    a_0 = 3.72,
    MM = 0.27333) "Pb(HS)2,aq";

  //1405
  constant DataRecord Pb_HS3_n1(
    R = Modelica.Constants.R/Pb_HS3_n1.MM,
    G_ref = -79375,
    H_ref = -120098,
    S_ref = 284.847,
    a1 = 5.1629e-05,
    a2 = 9351.99,
    a3 = -0.00012729,
    a4 = -154934,
    c1 = 166.0638,
    c2 = 280433,
    w_ref = 250200,
    z = -1,
    a_0 = 3.72,
    MM = 0.30639) "Pb(HS)3-";

  //1406
  constant DataRecord Mg_HSiO3_p1(
    R = Modelica.Constants.R/Mg_HSiO3_p1.MM,
    G_ref = -1477057,
    H_ref = -1613769,
    S_ref = -99.496,
    a1 = 2.6313e-06,
    a2 = -2611.99,
    a3 = 0.00034295,
    a4 = -105474,
    c1 = 153.9218,
    c2 = 195401,
    w_ref = 383970,
    z = 1,
    a_0 = 3.72,
    MM = 0.1014) "Mg(HSiO3)+";

  //1407
  constant DataRecord Ca_HSiO3_p1(
    R = Modelica.Constants.R/Ca_HSiO3_p1.MM,
    G_ref = -1574435,
    H_ref = -1686608,
    S_ref = -8.326,
    a1 = 4.4547e-06,
    a2 = -2166.77,
    a3 = 0.00032545,
    a4 = -107315,
    c1 = 128.8873,
    c2 = 153214,
    w_ref = 243970,
    z = 1,
    a_0 = 3.72,
    MM = 0.11717) "Ca(HSiO3)+";


end SolutesData;
