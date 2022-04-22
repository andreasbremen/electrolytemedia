within ElectrolyteMedia.Tests.LiquidPhase;
model partial_Titration "Example on how to use the Electrolyte Media library"
  /* This model illustrates how to use a medium from the Electrolyte Media library. It shows how
  to set a medium and how to implement corresponding mass balances.
  The model describes the titration of an analyte solution with a titrant solution.
  At t=0, the system consists of analyte solution only. A constant feed of 
  titrant solution is added over time.
  The implementation allows a generic model formulation for titration simulations. With the Medium 
  redeclaration and setting of initial values, a simulation may be perform, as provided in the 
  model NaOH_HCL_titration.
  You can use this model by drag and drop in a new model instance. Then, fill in the 
  user interface for the specification of the Medium and the initial values. The user interface 
  can be opened by double clicking on the model in the diagram view.*/

  /* 1. 
  First a medium is created. Since the model is generic, we instantiate the replaceable package 
  “Medium”. This medium package can later be redeclared as any liquid mixture medium. The 
  redeclarated medium then contains all parameters and for given species and considered 
  dissociation reactions.*/
  replaceable package Medium =
       ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid  annotation (choicesAllMatching = true);

  /* 2.
  Then we instantiate the model Medium.BaseProperties. This model calculates thermodynamic state
  variables (e.g. density). It is important to set correct starting conditions (Xredstart, pstart,
  and Tstart) to initialize the BaseProperties model. This will be discussed later.*/
  Medium.BaseProperties medium(Xredstart = Xred_system0, Tstart = T, pstart = p);

  /* 3.
  To calculate thermodynamic states we need two independent state variables 
  and independent mass fractions. In this case, the state variables are pressure and temperature, 
  which are set here. Any other combination of independent state variables is possible, too (e.g.,
  enthalpy and pressure). The independent mass fractions are discussed in step 5.*/
  parameter SI.Pressure p = 1e5;
  parameter SI.Temperature T = 298.15;

  /* 4.
  Next, the parameters for this model are defined. We provide the initial molalities, which will be 
  converted to mass fractions in step 5 (for mass balancing). Additionally, we define the initial 
  volume V_analyte. From the initial volume, we calculate the mass of the analyte m0_analyte. 
  Similarly, the titrant mass flow rate is defined.*/
  parameter Real[Medium.nL-1] molality_titrant;
  parameter Real[Medium.nL-1] molality_analyte;
  parameter SI.Volume V_analyte;
  final parameter SI.Mass m0_analyte = V_analyte*dl_analyte;
  final parameter SI.Density dl_analyte = Medium.Functions.calc_d(T,p,medium.Xfullstart);
  final parameter SI.MassFlowRate[Medium.nL] mdotfull_titrant = mDot_titrant*Xfull_titrant;
  parameter SI.MassFlowRate mDot_titrant;

  /* 5.
  There are two types of mass fraction vectors we use. The full mass fraction vector (Xfull) 
  includes the mass fraction of all species. The reduced mass fraction vector (Xred) contains 
  normalized independent reaction invariants. The origin of the reduced mass fraction vector is 
  further explained in detail in Bremen et al. (2021) and Kakhu and Pantelides (2003, 2004). For 
  now, it is relevant to know that the full mass fraction vector is used for the calculation of the 
  thermodynamic state variables and the reduced mass fraction vector. The reduced mass fraction 
  vector is used to set up mass balances.
  Naturally, we need to specify the full mass fraction vectors at t=0 without consideration of 
  equilibrium conditions, i.e., we provide the full mass fraction vector of non-dissociated 
  species. We then can calculate the initial reduced mass fraction vector of the system to pass to 
  the BaseProperties model in step 2. This will ensure that we use the correct initial masses for 
  the simulation. The BaseProperties model then considers the equilibrium conditions (isopotential,
  dissociation) and calculates the full mass fraction vector at equilibrium.*/
  final parameter SI.MassFraction[Medium.nL] Xfull_titrant = Medium.Functions.calc_Xfromm(molality_titrant);
  final parameter SI.MassFraction[Medium.nL] Xfull_analyte = Medium.Functions.calc_Xfromm(molality_analyte);
  final parameter Real[Medium.nX] Xred_system0= Medium.calc_Xred_from_Xfull(Xfull_analyte);
  final parameter SI.Mass[Medium.nL] mfull_system = m0_analyte*Xfull_analyte;

  /* 6.
  Next we create the time dependent variables we want to solve for. These are: the reduced mass 
  vector, the reduced mass fraction vector and the pH.*/
  SI.Mass[Medium.nX] mred_system;
  Real[Medium.nX] Xred_system;
  Real pH;

initial equation

  /*7.
  Here, weset the initial conditions for the system. In the model, we balacne the mass based 
  reaction invariants, and hence, we provide initial conditions for these.*/
  mred_system = Medium.calc_mred_from_mfull(mfull_system);

equation

  /* 8.
  The mass balance equations balance the mass based reaction invariants. With changing temperature 
  and pressure, the reaction invariants remain constant. A change of composition, however, 
  influences the reaction invariants (e.g., with a feed stream). The full feed mass flow vector is
  premultiplied by the transpose of the null space lambda_mass of the mass based stoichiometry 
  matrix nu_mass. Here, the reaction invariants mred_system changes by a constant feed of titrant.*/

  der(mred_system) = Medium.calc_mred_from_mfull(mdotfull_titrant);

  /* 9.
  We normalize the reduced mass fraction vector of the system. The reduced mass fraction vector is 
  required for the BaseProperties model.*/
  Xred_system = mred_system/sum(mred_system);

  /* 10.
  The base properties model needs the independent reduced mass fractions and two further 
  thermodynamic properties, e.g., pressure and temperature.*/
  medium.p = p;
  medium.T = T;
  medium.X[1:Medium.nX-1] = Xred_system[1:Medium.nX-1];

  /* 11.
  The pH is then calculated within the BaseProperties model. Further model variables are found 
  within the BaseProperties model.*/
  pH = medium.pH;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(StopTime=2000, __Dymola_NumberOfIntervals=4000),
    Documentation(info="<html>
</html>"));
end partial_Titration;
