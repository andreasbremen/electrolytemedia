within ElectrolyteMedia.Media.SolidPhase.Common;
partial package MixtureSolid "Medium model of a mixture of solid species based on Modelica library"
  extends Modelica.Icons.VariantsPackage;

  extends Modelica.Media.Interfaces.PartialMixtureMedium(
    ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX,
    redeclare final record FluidConstants =
        Modelica.Media.Interfaces.Types.IdealGas.FluidConstants,
    singleState=false,
    Temperature(min=273.15,max=573.15,start=300,nominal=300),
    Density(start=10, nominal=10),
    AbsolutePressure(start=10e5, nominal=10e5),
    mediumName=userInterface.mediumName,
    substanceNames = userInterface.substanceNames,
    T_default = userInterface.Tstart,
    p_default = userInterface.pstart,
    reference_X = userInterface.refX,
    nX = userInterface.nX);

  redeclare record extends ThermodynamicState "Thermodynamic state variables"
  end ThermodynamicState;
  constant DataRecord[nX] data = userInterface.data
    "Data record of solid species";

  import Modelica.Math;
  import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;

//   constant FluidConstants[nS] fluidConstants "Constant data for the fluid";

  constant MolarMass[nX] MMX= Functions.calc_MMX() "Vector of molar masses from each species";

  constant SI.Temperature T_0 = Common.T_0 "Temperature at reference point";
  constant SI.Pressure p_0 = Common.p_0 "Pressure at reference point";

  constant UserInterface userInterface  "User interface to create new medium" annotation (Placement(transformation(extent={{-10,-8},{10,12}})));

  redeclare model extends BaseProperties(
    T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default))
    "Base properties of ideal gas medium"
  equation
   assert(T >= 273.15 and T <= 573.15, "
  Temperature T (= " + String(T) + " K) is not in the allowed range
  0°C <= T <= 300°C required from medium model \"" + mediumName + "\".
  ");
    MM = molarMass(state);
    R = Functions.calc_R(X);
    h =Functions.calc_h(T,p,X);
    u =Functions.calc_u(T,p,X);

    // Has to be written in the form d=f(p,T) in order that static
    // state selection for p and T is possible
    d =Functions.calc_d(T,p,X);
    // connect state with BaseProperties
    state.T = T;
    state.p = p;
    state.X = X;
  end BaseProperties;

  redeclare function setState_pTX
    "Return thermodynamic state as function of p, T and composition X"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input Temperature T "Temperature";
      input MassFraction X[:]=reference_X "Mass fractions";
      output ThermodynamicState state;
  algorithm
      state := ThermodynamicState(p=p,T=T,X=X);
      annotation(Inline=true,smoothOrder=2);
  end setState_pTX;

  redeclare function setState_phX
    "Return thermodynamic state as function of p, h and composition X"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific enthalpy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
  algorithm
    state :=ThermodynamicState(p=p, T=T_hpX(h,p,X),X=X);
    annotation (Inline=true, smoothOrder=2);
  end setState_phX;

  redeclare function setState_psX
    "Return thermodynamic state as function of p, s and composition X"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific entropy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
  algorithm
    state :=ThermodynamicState(
              p=p,
              T=T_sX(s,X),
              X=X);
    annotation (Inline=true, smoothOrder=2);
  end setState_psX;

  redeclare function setState_dTX
    "Return thermodynamic state as function of d, T and composition X"
      extends Modelica.Icons.Function;
      input Density d "Density";
      input Temperature T "Temperature";
      input MassFraction X[:]=reference_X "Mass fractions";
      output ThermodynamicState state;
  algorithm
      state :=ThermodynamicState(
              p=p_dTX(
                T,
                d,
                X),
              T=T,
              X=X);
      annotation(Inline=true,smoothOrder=2);
  end setState_dTX;

  redeclare function extends setSmoothState
    "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
    //   algorithm
    //         state := ThermodynamicState(p=Modelica.Media.Common.smoothStep(x, state_a.p, state_b.p, x_small),
    //                                     T=Modelica.Media.Common.smoothStep(x, state_a.T, state_b.T, x_small));
    //         annotation(Inline=true,smoothOrder=2);
  end setSmoothState;

  redeclare function extends pressure "Return pressure of ideal gas"
  algorithm
    p := state.p;
    annotation(Inline=true,smoothOrder=2);
  end pressure;

  redeclare function extends temperature "Return temperature of ideal gas"
  algorithm
    T := state.T;
    annotation(Inline=true,smoothOrder=2);
  end temperature;

  redeclare function extends density "Return density of ideal gas"

  algorithm
    d := Functions.calc_d(state.T,state.p,state.X);
    annotation (Inline=true, smoothOrder=2);
  end density;

  redeclare function extends specificEnthalpy "Return specific enthalpy"
    extends Modelica.Icons.Function;

  algorithm
    h := Functions.calc_h(state.T,state.p,state.X);
  end specificEnthalpy;

  redeclare function extends specificInternalEnergy
    "Return specific internal energy"
    extends Modelica.Icons.Function;
  algorithm
    u := Functions.calc_u(state.T,state.p,state.X);

    annotation (Inline=true, smoothOrder=2);
  end specificInternalEnergy;

  redeclare function extends specificEntropy "Return specific entropy"
    extends Modelica.Icons.Function;
  algorithm
     s := Functions.calc_s(state.T,state.X);
    annotation (Inline=true, smoothOrder=2);
  end specificEntropy;

  redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
    extends Modelica.Icons.Function;
  algorithm
    g := Functions.calc_g(state.T,state.p,state.X);
    annotation (Inline=true, smoothOrder=2);
  end specificGibbsEnergy;

  redeclare function extends specificHelmholtzEnergy
    "Return specific Helmholtz energy"
    extends Modelica.Icons.Function;
  algorithm
    f := specificInternalEnergy(state) - state.T * specificEntropy(state);
    annotation(Inline=true,smoothOrder=2);
  end specificHelmholtzEnergy;

  redeclare function extends specificHeatCapacityCp
    "Return specific heat capacity at constant pressure"
  algorithm
    cp :=Functions.calc_cp(state.T,state.X);

    annotation(Inline=true,smoothOrder=2);
  end specificHeatCapacityCp;

  redeclare function extends specificHeatCapacityCv
    "Compute specific heat capacity at constant volume from temperature and gas data. cp equals cv in solid model"
  algorithm
    cv :=Functions.calc_cp(state.T,state.X);
    annotation(Inline=true,smoothOrder=2);
  end specificHeatCapacityCv;

  redeclare function extends molarMass "Return the molar mass of the medium"

  algorithm
    MM := Functions.calc_MM(state.X);

    annotation (Inline=true, smoothOrder=2);
  end molarMass;
end MixtureSolid;
