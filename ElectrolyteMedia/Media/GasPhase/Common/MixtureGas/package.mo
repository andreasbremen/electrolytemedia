within ElectrolyteMedia.Media.GasPhase.Common;
package MixtureGas "Medium model of a mixture of gases based on Modelica library"
  import Modelica.Math;
  import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;

  extends PartialMixtureMedium(
    redeclare replaceable record FluidConstants =
      Modelica.Media.Interfaces.Types.IdealGas.FluidConstants,
    ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.dTX,
    singleState = false,
    reducedX = false,
    reference_X=userInterface.refX,
    mediumName = userInterface.mediumName,
    substanceNames=userInterface.gasSubstanceNames,
    Density(start=0.10, nominal=0.10),
    AbsolutePressure(start=10e5, nominal=10e5),
    fluidConstants = userInterface.fluidconstants,
    Temperature(min=273.15, max=573.15, start=300, nominal=300),
     T_default = userInterface.Tstart,
     p_default = userInterface.pstart);

  redeclare record extends ThermodynamicState "Thermodynamic state variables"
  end ThermodynamicState;
  constant UserInterface userInterface annotation (Placement(transformation(extent={{-10,-8},{10,12}})));

  constant GasDataRecord[nX] data = userInterface.datag "Data record for single gas species parameters";
  constant InteractionDataRecord interaction = userInterface.interactiong "Data record containing coefficients for interactions beteween Peng-Robinson gases";

  constant Integer nG = nX "Number of gases";

  constant Media.Common.Types.GasModel GasModel = userInterface.GasModel "Enumeration of all gas types";

  constant MolarMass[nX] MMX= Functions.calc_MMX() "Returns vector containing molar mass of all species";

  redeclare function extends gasConstant "Return gasConstant"
  algorithm
    R :=Functions.calc_R(state.X);
    annotation(Inline = true, smoothOrder = 3);
  end gasConstant;

  redeclare replaceable model extends BaseProperties(
    T(stateSelect= StateSelect.always),
    d(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    Xi(each stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    final standardOrderComponents=true)
    "Base properties (p, d, T, h, u, R, MM, X, and Xi of NASA mixture gas"
                                       // if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    MoleFraction[nX] Y(start = reference_X./MMX/(sum(reference_X./MMX)));
    //     parameter Data dataset = data;
  equation
    //     assert(T >= 200 and T <= 6000, "
    // Temperature T (=" + String(T) + " K = 200 K) is not in the allowed range
    // 200 K <= T <= 6000 K
    // required from medium model \"" + mediumName + "\".");
    Y = massToMoleFractions(X,MMX);
    MM = molarMass(state);
    h =Functions.calc_h(
        T,
        d,
        X);
    R =Functions.calc_R(X);
    u = h - R*T;
    p =Functions.calc_p(
        T,
        d,
        X);
    // connect state with BaseProperties
    state.T = T;
    state.d = d;
    state.X = X;// if fixedX then reference_X else X;
  end BaseProperties;

  redeclare function setState_pTX
    "Return thermodynamic state as function of p, T and composition X"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
protected
    Density d;
  algorithm
    d :=d_TpX(
      T,
      p,
      X);
    state := ThermodynamicState(d=d,T=T,X=X);
    annotation (Inline=true, smoothOrder=2);
  end setState_pTX;

  redeclare function setState_phX
    "Return thermodynamic state as function of p, h and composition X"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEnthalpy h "Specific Enthalpy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
protected
    Temperature T;
    Density d;
  algorithm
    T :=T_phX(
      p,
      h,
      X);
    d :=d_TpX(
      T,
      p,
      X);
    state := ThermodynamicState(d=d,T=T,X=X);
    annotation (Inline=true, smoothOrder=2);
  end setState_phX;

  redeclare function setState_psX
    "Return thermodynamic state as function of p, s and composition X"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input SpecificEntropy s "Specific Enthalpy";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
protected
    Temperature T;
    Density d;
  algorithm
    T :=T_psX(
      p,
      s,
      X);
    d :=d_TpX(
      T,
      p,
      X);
    state := ThermodynamicState(d=d,T=T,X=X);
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
    state := ThermodynamicState(d=d, T=T, X=X);
    annotation (Inline=true, smoothOrder=2);
  end setState_dTX;

  redeclare function extends pressure "Return pressure of PR gas"
  algorithm
    p :=Functions.calc_p(
        state.T,
        state.d,
        state.X);
    annotation (Inline=true, smoothOrder=2);
  end pressure;

  redeclare function extends temperature "Return temperature of ideal gas"
  algorithm
    T := state.T;
    annotation(Inline=true,smoothOrder=2);
  end temperature;

  redeclare function extends density "Return density of ideal gas"
  algorithm
    d := state.d;
    annotation(Inline = true, smoothOrder = 3);
  end density;

  redeclare function extends molarMass "Return molar mass of mixture"
  algorithm
    MM :=Functions.calc_MM(state.X);
    annotation(Inline=true,smoothOrder=2);
  end molarMass;

  redeclare function extends specificEnthalpy "Return specific enthalpy"
    extends Modelica.Icons.Function;
  algorithm
    h :=Functions.calc_h(
        state.T,
        state.d,
        state.X);
    annotation(Inline=true,smoothOrder=2);
  end specificEnthalpy;

  redeclare function extends specificInternalEnergy "Return specific internal energy"
    extends Modelica.Icons.Function;
  algorithm
    u :=Functions.calc_u(
        state.T,
        state.d,
        state.X);
    annotation(Inline=true,smoothOrder=2);
  end specificInternalEnergy;

  redeclare function extends specificEntropy "Return specific entropy"
    extends Modelica.Icons.Function;
  algorithm
    s :=Functions.calc_s(
        state.T,
        state.d,
        state.X);
    annotation(Inline=true,smoothOrder=2);
  end specificEntropy;

  redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
    extends Modelica.Icons.Function;
  algorithm
    g :=Functions.calc_g(
        state.T,
        state.d,
        state.X);
    annotation(Inline=true,smoothOrder=2);
  end specificGibbsEnergy;

  redeclare function extends specificHelmholtzEnergy
    "Return specific Helmholtz energy"
    extends Modelica.Icons.Function;
  algorithm
    f := specificEnthalpy(state) - gasConstant(state)*state.T - state.T*specificEntropy(state);
    annotation(Inline=true,smoothOrder=2);
  end specificHelmholtzEnergy;

  redeclare function extends specificHeatCapacityCp
    "Return specific heat capacity at constant pressure"
  algorithm
    cp :=Functions.calc_cp(
        state.T,
        state.d,
        state.X);
    annotation(Inline=true,smoothOrder=1);
  end specificHeatCapacityCp;

  redeclare function extends specificHeatCapacityCv
    "Compute specific heat capacity at constant volume from temperature and gas data"
  algorithm
    cv :=Functions.calc_cv(
        state.T,
        state.d,
        state.X);
    annotation(Inline=true,smoothOrder=2);
  end specificHeatCapacityCv;

  redeclare function extends isentropicExponent "Return isentropic exponent"
  algorithm
    gamma := specificHeatCapacityCp(state)/specificHeatCapacityCv(state);
    annotation(Inline=true,smoothOrder=2);
  end isentropicExponent;

  redeclare function extends velocityOfSound "Return velocity of sound"
    extends Modelica.Icons.Function;
protected
    MolarVolume v = molarVolume(state);
    MolarMass MM = molarMass(state);
    Real kappa(unit="1/Pa") = isothermalCompressibility(state);
  algorithm
    a := sqrt(max(0,v/(kappa*MM)));
    annotation(Inline=true,smoothOrder=2);
  end velocityOfSound;

  redeclare function extends isobaricExpansionCoefficient
    "Returns overall the isobaric expansion coefficient beta"
  algorithm
    beta :=Functions.calc_beta(
        state.T,
        state.d,
        state.X);
    annotation(Inline=true,smoothOrder=2);
  end isobaricExpansionCoefficient;

  redeclare function extends isothermalCompressibility
    "Returns overall the isothermal compressibility factor"
  algorithm
    kappa :=Functions.calc_kappa(
        state.T,
        state.d,
        state.X);
    annotation(Inline=true,smoothOrder=2);
  end isothermalCompressibility;

  redeclare function extends density_derp_T
    "Returns the partial derivative of density with respect to pressure at constant temperature"
protected
    MolarMass MM = molarMass(state);
  algorithm
  ddpT := ElectrolyteMedia.Media.GasPhase.Common.Functions.calc_der_d_p(
      state.T,
      state.d,
      state.X);
    annotation(Inline=true,smoothOrder=2);
  end density_derp_T;

  redeclare function extends density_derT_p
    "Returns the partial derivative of density with respect to temperature at constant pressure"
protected
    MolarMass MM = molarMass(state);
  algorithm
  ddTp := ElectrolyteMedia.Media.GasPhase.Common.Functions.calc_der_d_T(
      state.T,
      state.d,
      state.X);
    annotation(Inline=true,smoothOrder=2);
  end density_derT_p;

  redeclare replaceable function extends dynamicViscosity
    "Return mixture dynamic viscosity"
protected
  DynamicViscosity[nX] etaX "Component dynamic viscosities";
  algorithm
  for i in 1:nX loop
      etaX[i] := dynamicViscosityLowPressure(        state.T,
                   fluidConstants[i].criticalTemperature,
                   fluidConstants[i].molarMass,
                   fluidConstants[i].criticalMolarVolume,
                   fluidConstants[i].acentricFactor,
                   fluidConstants[i].dipoleMoment);
  end for;
  eta := gasMixtureViscosity(massToMoleFractions(state.X,
                         fluidConstants[:].molarMass),
             fluidConstants[:].molarMass,
             etaX);
  annotation (smoothOrder=2);
  end dynamicViscosity;
end MixtureGas;
