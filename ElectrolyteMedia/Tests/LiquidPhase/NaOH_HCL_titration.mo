within ElectrolyteMedia.Tests.LiquidPhase;
model NaOH_HCL_titration "Implementation of partial_Titration"

  /*This model describes titration of a Sodium Hydroxide solution (NaOH solution as Na+ and 
  OH- ions) with a Hydrogen Chloride (HCl) solution. Sodium Hydroxide is the analyte and 
  Hydrogen Chloride the titrant. At t=0, the system conists of analyte (NaOH solution) only. 
  A constant feed of titrant (HCl solution) is added to the analyte.
  To set up a simulation, all parameters have to be provided, i.e., the Medium needs to be 
  redeclared and the further model parameters need to be set.
  
  First, pick a medium that describes the two solutions, i.e., set a medium that includes both
  solutions as a mixture.
  Present species:        HCl, Na+, Cl-, OH-, H+, H2O
  Considered reactions:   HCl <-> Cl- + H+
                          H2O <-> OH- + H+ 
  It is assumed that NaOH is completely dissociated. Thus, there is no equilibrium reaction and 
  non-dissociated NaOH.
  A medium that meets these requirements is the "HydrochloricAcid" medium. Hence, we redeclare the 
  Medium package with the HydrochloricAcid medium.
  
  Secondly, the initial molalities of the analyte and the constant molalities of the titrant are 
  provided. Lastly, the initial volume of the analyte and the mass flow rate of the titrant is
  specified.*/

  partial_Titration partial_Titration1(
    redeclare package Medium =
        Media.LiquidPhase.MixtureLiquids.HydrochloricAcid,
    molality_titrant={0,0,0.35,0,0.35},
    molality_analyte={0,0.24,0,0.24,0},
    V_analyte=25e-6,
    mDot_titrant=2.5e-5) annotation (Placement(transformation(extent={{-10,-14},{10,6}})));
   annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=2000, __Dymola_NumberOfIntervals=4000));
end NaOH_HCL_titration;
