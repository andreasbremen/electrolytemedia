within ElectrolyteMedia.Media.SolidLiquidPhase.Common;
package MixtureLiquid
  import Modelica.Math;
  import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;

  extends Modelica.Media.Interfaces.PartialMixtureMedium(
     ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX,
     reference_X=userInterface.refX,
     nX = nS-nR,
     reducedX = true,
     nXi = nX-1,
     mediumName = userInterface.mediumName,
     singleState=false,
     substanceNames=cat(1,userInterface.solidSubstanceNames,userInterface.liquidSubstanceNames),
     Density(start=10, nominal=10),
     AbsolutePressure(start=10e5, nominal=10e5),
     Temperature(min=273.15, max=573.15, start=300, nominal=300),
     T_default = userInterface.Tstart,
     p_default = userInterface.pstart);

  redeclare record extends ThermodynamicState "Thermodynamic state variables"
    Real[nF] Xfull;
  end ThermodynamicState;

  constant SI.Pressure prefg = 1e5 "Reference pressure for gas phase";
  constant DataRecordL[nL-1] datal=userInterface.datal;
  constant LiquidInteractionDataRecord interactionL=userInterface.interactionL;
  constant DataRecordS[ns] datas=userInterface.datas;

  constant Integer nL(min=2)=userInterface.nL;
  constant Integer nLi(min=1) = nL-1;
  constant Integer ns(min=1)=userInterface.ns;
  constant Integer nR=userInterface.nR;
  constant Integer nF = nX + nR;

  constant Media.Common.Types.LiquidModel LiquidModel=userInterface.LiquidModel;
  constant Real[nR,nF] nu = userInterface.nu;
  constant Real[nF,nX] lambda= Media.Common.Reaction.calc_lambda(nu, nR, nF) "Nullspace of stoichiometry matrix";
  constant Real[nF,nX] lambda_mass={{lambda[i,j]/MMX[i] for j in 1:nX} for i in 1:nF} "Nullspace of mass based stoichiometry matrix";

  constant Real[nR,nF] nu_id= Media.Common.Reaction.calc_nu_id(nu) "Identity transformation of stoichiometry matrix for initialization";
  constant Integer[nF,nF] P_to_orig=Media.Common.Reaction.calc_P_nu_id(nu) "Permutation matrix between original order and identity transformation";
  constant Integer[nF,nF] P_to_id = transpose(P_to_orig) "Permutation matrix between identity transformation and original order";
  constant Real[nF,nX] lambda_id = Media.Common.Reaction.calc_lambda_id(nu_id) "lambda in identity transformation";
                                                                                                                   // transpose(cat(2,-transpose(nu_id[:,nR+1:nF]),identity(nX)))
  constant Real[nF] MMX_id = P_to_id*MMX;

  constant UserInterface
    userInterface "User interface"
    annotation (Placement(transformation(extent={{-10,-8},{10,12}})));


  constant MolarMass MH2O = IF97.MH2O "Molar mass of solvent";
  constant SpecificHeatCapacity RH2O = IF97.RH2O;

  constant MolarMass[ns + nL] MMX=Functions.calc_MMX();

  constant Integer Hindex = max(2,Modelica.Math.BooleanVectors.firstTrueIndex({abs(datal[i].MM-0.001008) == 0 for i in 1:nL-1})+ns);
  constant Integer OHindex = Modelica.Math.BooleanVectors.firstTrueIndex({abs(datal[i].MM-0.017008) == 0 for i in 1:nL-1})+ns;

  redeclare function extends gasConstant "Return gasConstant"
  algorithm
    R :=Functions.calc_R(state.Xfull);
    annotation(Inline = true, smoothOrder = 3);
  end gasConstant;

  redeclare replaceable model extends BaseProperties(
    T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    Xi(each stateSelect=if preferredMediumStates then StateSelect.always else StateSelect.default),
    final standardOrderComponents = true)
    "Base properties (p, d, T, h, u, R, MM, X, and Xfull of aqueous liquid solution"
    MoleFraction[ns+nL] Yfull(start=Yfullstart) "Full mole fraction vector";
    MassFraction[ns+nL] Xfull(start=Xfullstart,each min = 0) "Full mass fraction vector";
    SI.Mass[ns+nL] mfull(start=mfullstart*scale)  "Full amount of substance vector (substitute)";
    SI.SpecificEnergy[ns+nL] gi;
    SI.MolarEnergy[ns+nL] gim;
    SI.MolarEnergy[nR] gr "Gibbs free energy of reaction at T and p ";
    Real[nR] K "Equilibrium constant";
    Real[nR] logK "Logarithmic equilibrium constant";
    Real[nR] log10K "Decadic logarithmic equilibrium constant";
    Real[nL] m "Molality";
    Real[nL] a "Activity";
    Real[nL] gamma "Activity coefficient";
    Real[nL] loga( start = logastart) "Logarithmic activity of all species";
    Real[ns+nL] logreacBase(start = cat(1,-zstart[1:ns],logastart));//-zstart[ns+1:ns+nL]));
    Real pH;
    Real I;
    Real z [ns](start = zstart[1:ns]);
    MoleFraction[ns] Ys;
    MoleFraction[nL] Yl;
    MassFraction[ns] Xs;
    MassFraction[nL] Xl;
    SI.VolumeFraction Phil;
    SI.Density dl;

    parameter MoleFraction[ns + nL] Yfullstart=Functions.calc_Y_M(Xfullstart, MMX);
    parameter Real[ns+nL] zstart = {1e-20./(2*max(1e-20,Xfullstart[i])) for i in 1:ns+nL};

  parameter Real[ns + nL] Xfullstart=X_pTXred(
        pstart,
        Tstart,
        Xredstart);
    parameter MassFraction[nX] Xredstart = reference_X;
    parameter Real[nX] Xredstart_ = transpose(lambda_mass)*Xfullstart "for scaling of species moles consistent with maxx invariant vector X";
    parameter Real[ns+nL] mfullstart = Xfullstart/sum(Xredstart_);
    final parameter Real[nL-1] mstart_ = Functions.LiquidFunctions.calc_mfromX(Xlstart);
    final parameter Real[nL] mstart = cat(1,mstart_,{Ylstart[nL]/MH2O});
    final parameter Real[nL] gammastart = Functions.LiquidFunctions.calc_gamma(Tstart,pstart,Xlstart);
    final parameter Real[nL] astart = gammastart.*mstart;
    final parameter Real[nL] logastart = log(astart);
    final parameter Real[ns+nL] logastart_ = cat(1,zeros(ns),logastart);
    parameter MassFraction[ns] Xsstart = Xfullstart[1:ns]/sum(Xfullstart[1:ns]);
    parameter MassFraction[nL] Xlstart = Xfullstart[1+ns:ns+nL]/sum(Xfullstart[1+ns:ns+nL]);
    final parameter MoleFraction[nL] Ylstart = Functions.LiquidFunctions.calc_Y(Xlstart);
    parameter Temperature Tstart = T_default;
    parameter SI.Pressure pstart = p_default;
    parameter Real eps= 1e-10;

    parameter Real[:,:] lambda_mass_scaled = lambda_mass/scale;
    parameter Real scale = 100000;//1000000000;//

    Real[nX] X_ "Substitute for X to ensure non-reactive species to be > 0";
  equation
    assert(T >= 273.15 and T <= 573.15, "Temperature T (=" + String(T) + " K = 200 K) is not in the allowed range 273.15 K <= T <= 573.15 K required from medium model \"" + mediumName + "\".");
    Yfull = massToMoleFractions(Xfull,MMX);
    Xs[1:ns-1] = Xfull[1:ns-1]/sum(Xfull[1:ns]);
    sum(Xs) = 1;
    Xl[1:nL-1] = Xfull[1+ns:ns+nL-1]/sum(Xfull[1+ns:ns+nL]);
    sum(Xl) = 1;
    Ys = massToMoleFractions(Xs,MMX[1:ns]);
    Yl = massToMoleFractions(Xl,MMX[1+ns:ns+nL]);
    MM = molarMass(state);
    h =Functions.calc_h(T,p,Xfull);
    R =Functions.calc_R(Xfull);
    u =Functions.calc_u(T,p,Xfull);
    d =Functions.calc_d(T,p,Xfull);
    // connect state with BaseProperties
    state.T = T;
    state.p = p;
    state.X = X;
    state.Xfull = Xfull;

    Xfull[1:nF] = mfull[1:nF]/sum(mfull);
    X_ = transpose(lambda_mass_scaled)*mfull;
    X_ = {max(1e-25,X[i]) for i in 1:nX};

    gi =Functions.calc_gi(T, p);
    gim = gi.*MMX;
    m[1:nL-1] =Functions.LiquidFunctions.calc_mfromY(Yl);
    m[nL] = Yl[nL]/IF97.MH2O;
    gamma = Functions.LiquidFunctions.calc_gamma(T,p,Xl);
    a = gamma.*m;
    I =Functions.LiquidFunctions.calc_I(Xl);
    gr = nu*gim;
    logK = -gr / Modelica.Constants.R / T;
    K = exp(logK);
    log10K = log10(K);
    logreacBase = cat(1,-z[1:ns],loga);//-z[1+ns:ns+nL]);
    //    nu*logreacBase = logK;
     nu*logreacBase./logK = ones(nR);
    a = exp(loga);
    pH = -loga[Hindex-ns]*log10(exp(1));
    Phil = Functions.calc_Phil(T,p,Xfull);
    dl = Functions.LiquidFunctions.calc_d(T,p,Xl);

    //Fischer-Burmeister function
    zeros(ns) = mfull[1:ns] + z - sqrt(mfull[1:ns].^2+z.^2 .+ eps);

    annotation (experiment(
        StopTime=3000,
        __Dymola_NumberOfIntervals=5000,
        Tolerance=1e-10,
        __Dymola_Algorithm="Euler"));
  end BaseProperties;

  redeclare function setState_pTX
    "Return thermodynamic state as function of p, T and composition X"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:] "Mass fractions";
    output ThermodynamicState state;
protected
    MassFraction[nX] Xred = if size(X,1)==nX then X else cat(1,X,{1-sum(X)});
  Real[nF] init=X_pTXred(p,T,Xred);
    MassFraction[nF] Xfull = init[1:nF];
  algorithm
    state := ThermodynamicState(p=p,T=T,X=Xred,Xfull=Xfull);//,Yred=Yred);
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
    MassFraction[nX] Xred = if size(X,1)==nX then X else cat(1,X,{1-sum(X)});
    MassFraction[nF] Xfull;
    Temperature T;
  algorithm
       (Xfull,T) :=XT_phXred(
         p,
         h,
         Xred);
    state := ThermodynamicState(p=p,T=T,X=Xred,Xfull=Xfull);
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
    MassFraction[nX] Xred = if size(X,1)==nX then X else cat(1,X,{1-sum(X)});
    MassFraction[nF] Xfull;
    Temperature T;
    MassFraction[ns] Xg;
    SI.Density dg "Gas density";
  algorithm
    (Xfull,T) :=XT_psXred(p,s,Xred);
    state := ThermodynamicState(p=p,T=T,X=Xred,Xfull = Xfull);
    annotation (Inline=true, smoothOrder=2);
  end setState_psX;

  redeclare function setState_dTX
    "Return thermodynamic state as function of d, T and composition X"
    extends Modelica.Icons.Function;
    input Density d "Density";
    input Temperature T "Temperature";
    input MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
protected
    MassFraction[nX] Xred = if size(X,1)==nX then X else cat(1,X,{1-sum(X)});
    MassFraction[nF] Xfull;
    SI.Pressure p;
  algorithm
    (Xfull,p) :=Xp_dTXred(d,T,Xred);
    state := ThermodynamicState(p=p,T=T,X=Xred,Xfull = Xfull);
    annotation (Inline=true, smoothOrder=2);
  end setState_dTX;

  redeclare function extends pressure "Return pressure of PR gas"
  algorithm
    p := state.p;
    annotation (Inline=true, smoothOrder=2);
  end pressure;

  redeclare function extends temperature "Return temperature of ideal gas"
  algorithm
    T := state.T;
    annotation(Inline=true,smoothOrder=2);
  end temperature;

  redeclare function extends density "Return density of ideal gas"
  algorithm
    d :=Functions.calc_d(
      state.T,
      state.p,
      state.Xfull);
    annotation(Inline = true, smoothOrder = 3);
  end density;


  redeclare function extends molarMass "Return molar mass of mixture"
  algorithm
    MM := Functions.calc_MM(state.Xfull);
    annotation(Inline=true,smoothOrder=2);
  end molarMass;

  redeclare function extends specificEnthalpy "Return specific enthalpy"
    extends Modelica.Icons.Function;
  algorithm
    h :=Functions.calc_h(
      state.T,
      state.p,
      state.Xfull);
    annotation(Inline=true,smoothOrder=2);
  end specificEnthalpy;

  redeclare function extends specificInternalEnergy "Return specific internal energy"
    extends Modelica.Icons.Function;
  algorithm
    u :=Functions.calc_u(
      state.T,
      state.p,
      state.Xfull);
    annotation(Inline=true,smoothOrder=2);
  end specificInternalEnergy;

  redeclare function extends specificEntropy "Return specific entropy"
    extends Modelica.Icons.Function;
  algorithm
    s :=Functions.calc_s(
      state.T,
      state.p,
      state.Xfull);
    annotation(Inline=true,smoothOrder=2);
  end specificEntropy;

  redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
    extends Modelica.Icons.Function;
  algorithm
    g :=Functions.calc_g(
      state.T,
      state.p,
      state.Xfull);
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
      state.p,
      state.Xfull);
    annotation(Inline=true,smoothOrder=1);
  end specificHeatCapacityCp;

  redeclare function extends specificHeatCapacityCv
    "Compute specific heat capacity at constant volume from temperature and gas data"
  algorithm
    cv :=Functions.calc_cv(
      state.T,
      state.p,
      state.Xfull);
    annotation(Inline=true,smoothOrder=2);
  end specificHeatCapacityCv;

  redeclare function extends isentropicExponent "Return isentropic exponent"
  algorithm
    gamma := specificHeatCapacityCp(state)/specificHeatCapacityCv(state);
    annotation(Inline=true,smoothOrder=2);
  end isentropicExponent;

  redeclare replaceable function extends dynamicViscosity
    "Return mixture dynamic viscosity"
  algorithm
  eta := 1;//Common.Functions.IF97_R1_Tp.calc_eta(state.p,state.T);
  annotation (smoothOrder=2);
  end dynamicViscosity;












end MixtureLiquid;
