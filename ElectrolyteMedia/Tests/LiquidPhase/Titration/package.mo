within ElectrolyteMedia.Tests.LiquidPhase;
package Titration
  model TitrationModel
   replaceable package Medium =
        ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid
        annotation (choicesAllMatching = true);

    Medium.BaseProperties medium_titrant(Xredstart = Xred_titrant);
    Medium.BaseProperties medium_solution(Xredstart = Xred_analyte);

    final parameter SI.Pressure p = 1e5;
    final parameter SI.Temperature T = 298.15;

    parameter Real[Medium.nL-1] molality_titrant;
    parameter Real[Medium.nL-1] molality_analyte;
    parameter SI.Volume V_analyte;

    final parameter SI.MassFraction[Medium.nL] Xfull_titrant = Medium.Functions.calc_Xfromm(molality_titrant);
    final parameter SI.MassFraction[Medium.nL] Xfull_analyte = Medium.Functions.calc_Xfromm(molality_analyte);
    final parameter Real[Medium.nX] Xred_titrant = transpose(Medium.lambda_mass)*Xfull_titrant/(sum(transpose(Medium.lambda_mass)*Xfull_titrant));
    final parameter Real[Medium.nX] Xred_analyte = transpose(Medium.lambda_mass)*Xfull_analyte/(sum(transpose(Medium.lambda_mass)*Xfull_analyte));

    SI.Volume V_titrant;
    SI.Mass m_titrant;
    SI.Mass[Medium.nX] mred;
    Real[Medium.nX] Xred_solution;

    Real pH;

    final parameter SI.MassFlowRate mDot = 2.5e-3;
    final parameter SI.Mass m0_analyte = V_analyte*dl_analyte;
    final parameter SI.Density dl_titrant = Medium.Functions.calc_d(T,p,medium_titrant.Xfullstart);
    final parameter SI.Density dl_analyte = Medium.Functions.calc_d(T,p,medium_solution.Xfullstart);

  initial equation

    m_titrant = 0;
    mred = m0_analyte*transpose(Medium.lambda_mass)*Xfull_analyte;

  equation

    der(m_titrant) = mDot;
    der(mred) = mDot*transpose(Medium.lambda_mass)*Xfull_titrant;

    m_titrant = V_titrant * dl_titrant;

    medium_solution.pH = pH;
    medium_solution.X[1:Medium.nX-1] = Xred_solution[1:Medium.nX-1];
    medium_solution.T = T;
    medium_solution.p = p;
    medium_titrant.T = T;
    medium_titrant.p = p;
    medium_titrant.X[1:Medium.nX-1] = Xred_titrant[1:Medium.nX-1];

    Xred_solution = mred/sum(mred);

    annotation (experiment(
        StopTime=1000,
        __Dymola_NumberOfIntervals=5000,
        Tolerance=1e-14,
        __Dymola_Algorithm="Dassl"));
  end TitrationModel;

  model SimulateTitration_HCl_NaOH
    TitrationModel titration_Model(
      redeclare package Medium =
          Media.LiquidPhase.MixtureLiquids.HydrochloricAcid,
      molality_titrant={0,0,0.34,0,0.34},
      molality_analyte={0,0.25,0,0.25,0},
      V_analyte=25e-6)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=2000,
        __Dymola_NumberOfIntervals=4000,
        __Dymola_Algorithm="Dassl"));
  end SimulateTitration_HCl_NaOH;

  model SimulateTitration_PhosphoricAcid_NaOH
    TitrationModel titration_Model(
      redeclare package Medium =
          Media.LiquidPhase.MixtureLiquids.PhosphoricAcid_NaOH,
      molality_titrant={0,0,0,0,0.1,0.1,0},
      molality_analyte={0.2,0,0,0,0,0,0},
      V_analyte=25e-6)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=60,
        __Dymola_NumberOfIntervals=1000,
        __Dymola_Algorithm="Dassl"));
  end SimulateTitration_PhosphoricAcid_NaOH;

  model SimulateTitration_OxcalicAcid
    TitrationModel titration_Model(
      redeclare package Medium =
          Media.LiquidPhase.MixtureLiquids.OxalicAcid_NaOH,
      molality_titrant={0,0,0,0.1,0.1,0},
      molality_analyte={0.2,0,0,0,0,0},
      V_analyte=25e-6)
      annotation (Placement(transformation(extent={{-10,-12},{10,8}})));

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(StopTime=10000, Tolerance=1e-08));
  end SimulateTitration_OxcalicAcid;

  model SimulateTitration_AceticAcid_NaOH
    TitrationModel titration_Model(
      redeclare package Medium =
          Media.LiquidPhase.MixtureLiquids.AceticAcid_NaOH,
      molality_titrant={0,0,0.15,0,0.15,0},
      molality_analyte={0.1,0,0,0,0,0},
      V_analyte=25e-6)
      annotation (Placement(transformation(extent={{-10,-12},{10,8}})));

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=1000,
        Tolerance=1e-08,
        __Dymola_Algorithm="Dassl"));
  end SimulateTitration_AceticAcid_NaOH;

  model SimulateTitration_AceticAcid_NaOH_dis
    TitrationModel titration_Model(
      redeclare package Medium =
          Media.LiquidPhase.MixtureLiquids.AceticAcid_NaOH_dis,
      molality_titrant={0,0,0.15,0,0,0,0},
      molality_analyte={0.1,0,0,0,0,0,0},
      V_analyte=25e-6)
      annotation (Placement(transformation(extent={{-10,-12},{10,8}})));

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false)),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=1000,
        Tolerance=1e-08,
        __Dymola_Algorithm="Dassl"));
  end SimulateTitration_AceticAcid_NaOH_dis;
end Titration;
