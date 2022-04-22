within ElectrolyteMedia.Media.GasLiquidPhase.Common;
package SaturatedLiquid "Mixture liquid at bubble point temperature with only H2O in gas phase. Temperature and pressure are not independent, and vapor fraction has to be provided"
  import Modelica.Math;
  import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;

  extends Modelica.Media.Interfaces.PartialMixtureMedium(
     ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX,
     fixedX=false,
     reducedX = true,
     singleState=false,
     reference_X=userInterface.refX,
     nX = nS-nR,
     nXi = nX-1,
     mediumName = userInterface.mediumName,
     substanceNames=cat(1,userInterface.gasSubstanceNames,userInterface.liquidSubstanceNames),
     Density(start=10, nominal=10),
     AbsolutePressure(start=10e5, nominal=10e5),
     Temperature(min=273.15, max=573.15, start=300, nominal=300),
     T_default = userInterface.Tstart,
     p_default = userInterface.pstart);

  redeclare record extends ThermodynamicState "Thermodynamic state variables"
    Real[nF] Xfull;
    SI.Density dg "Gas density";
  end ThermodynamicState;
  constant SI.Pressure prefg = 1e5 "Reference pressure for gas phase";
  constant DataRecordL[nL-1] datal=userInterface.datal;
  constant DataRecordG[nG] datag=userInterface.datag;

  constant GasInteractionDataRecord interactionG=userInterface.interactionG;
  constant LiquidInteractionDataRecord interactionL=userInterface.interactionL;

  constant Integer nL(min=2)=userInterface.nL;
  constant Integer nLi(min=1) = nL-1;
  constant Integer nG(min=1)=userInterface.nG;
  constant Integer nR=userInterface.nR;
  constant Integer nF = nX + nR;
  constant Integer nR_L = userInterface.nR_L;
  constant Integer nX_L = nL - nR_L;

  constant Media.Common.Types.LiquidModel LiquidModel=userInterface.LiquidModel;
  constant Media.Common.Types.GasModel GasModel=userInterface.GasModel;
  constant Real[nR,nF] nu = userInterface.nu;
  constant Real[nR,nF] nu_mass = Functions.Reaction.calc_nu_mass(nu);
  constant Real[nR_L,nL] nu_L = userInterface.nu_L;
  constant Real[nR_L,nL] nu_mass_L = Functions.Reaction.calc_nu_mass_L(nu_L);

  constant UserInterface userInterface  annotation (Placement(transformation(extent={{-10,-8},{10,12}})));


constant Real[:,:] lambda=
    Media.Common.Reaction.calc_lambda(
    nu,
    nR,
    nF);
constant Real[nF,nX] lambda_mass={{lambda[i,j]/MMX[i] for j in 1:nX} for i in 1:nF}*1000 "Null space of mass based stoichiometry matrix";

  constant Real[:,:] lambda_L = Media.Common.Reaction.calc_lambda(nu_L,nR_L,nL);
  constant Real[nL,nX_L] lambda_mass_L = {{lambda_L[i,j]/MMX[i+nG] for j in 1:nX_L} for i in 1:nL}*1000 "Null space of mass based stoichiometry matrix";
                                                                                                                                                        //Media.Common.Reaction.calc_lambda(nu_mass_L,nR_L,nL);//

  constant MolarMass MH2O = IF97.MH2O "Molar mass of solvent";
  constant SpecificHeatCapacity RH2O = IF97.RH2O;

constant MolarMass[nG + nL] MMX= Functions.calc_MMX() "Molar masses of components";

  constant Integer Hindex = max(2,Modelica.Math.BooleanVectors.firstTrueIndex({abs(datal[i].MM-0.001008) == 0 for i in 1:nL-1})+nG);
  constant Integer OHindex = Modelica.Math.BooleanVectors.firstTrueIndex({abs(datal[i].MM-0.017008) == 0 for i in 1:nL-1})+nG;

  redeclare function extends gasConstant "Return gasConstant"
  algorithm
    R :=Functions.calc_R(state.Xfull);
    annotation(Inline = true, smoothOrder = 3);
  end gasConstant;

  redeclare replaceable model extends BaseProperties(
    T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default, start = Tinitial),
    p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default, start = pstart),
    Xi(each stateSelect=if preferredMediumStates then StateSelect.always else StateSelect.default),
    final standardOrderComponents=true)
    "Base properties (p, d, T, h, u, R, MM, X, and Xfull of aqueous liquid solution"
    MoleFraction[nG+nL] Yfull(start=Yfullstart) "Full mole fraction vector";
    MassFraction[nG+nL] Xfull(start=Xfullstart,each min = 0) "Full mass fraction vector";
    SI.Mass[nG+nL] mfull(start=mfullstart)  "Full amount of substance vector (substitute)";
    SI.MolarEnergy[nG+nL] gim;
    SI.SpecificEnergy[nG+nL] gi;
    SI.MolarEnergy[nR] gr "Gibbs free energy of reaction at T and p ";
    Real[nR] K "Equilibrium constant";
    Real[nR] logK "Logarithmic equilibrium constant";
    Real[nL] m "Molality";
    Real[nL] a "Activity";
    Real[nL] loga( start = logastart);
    Real[nL] gamma "Activity coefficient";
    Real[nG+nL] logreacBase(start= cat(1,log(fstart),log(mstart))) "Logarithmic fugacity for gas and activity for liquid species";
    Real pH;
    Real I;
    MoleFraction[nG] Yg;
    MoleFraction[nL] Yl;
    MassFraction[nG] Xg;
    MassFraction[nL] Xl;
    SI.Density dg(start=dgstart) "Gas density";
    Real[nG] phi(start = ones(nG)) "Fugacity coefficient";
    Real[nG] f( start=fstart) "Fugacity";
    Real[nG] logf( start=log(fstart)) "Logarithmic Fugacity";
    Real vfrac(start = vfracstart) "Mass based vapor fraction";
    final parameter MoleFraction[nG + nL] Yfullstart = Functions.calc_Y_M(Xfullstart, MMX);  // Common.Functions.Reaction.solve_LE(mstart_,Tstart,pstart,data,nu,nF,nR);
    parameter MassFraction[nG + nL] Xfullstart=X_vfracpXred(
                vfracstart,
                pstart,
                Xredstart);                                                                                                                                //cat(1,fill(1/nG,nG),fill(1e-3,nL-1),{1-(nL-1)*1e-3});//
    parameter MassFraction[nX] Xredstart = reference_X;
    final parameter Real[nX] Xredstart_ = transpose(lambda_mass)*Xfullstart "for scaling of species moles consistent with maxx invariant vector X";
    final parameter Real[nG+nL] mfullstart = Xfullstart/sum(Xredstart_);
    final parameter Real[nL-1] mstart_ = Functions.LiquidFunctions.calc_mfromX(Xlstart);
    final parameter Real[nL] mstart = cat(1,mstart_,{Ylstart[nL]/MH2O});
    final parameter Real[nL] gammastart = Functions.LiquidFunctions.calc_gamma(Tstart,pstart,Xlstart);
    final parameter Real[nL] astart = gammastart.*mstart;
    final parameter Real[nL] logastart = log(astart);
    final parameter MassFraction[nG] Xgstart = Xfullstart[1:nG]/sum(Xfullstart[1:nG]);
    final parameter MassFraction[nL] Xlstart = Xfullstart[1+nG:nG+nL]/sum(Xfullstart[1+nG:nG+nL]);
    final parameter MoleFraction[nL] Ylstart = Functions.LiquidFunctions.calc_Y(Xlstart);
    final parameter SI.Density dgstart = dg_TpX(Tstart,pstart,Xgstart) "Start value for gas density";
    final parameter Real[nG] fstart = Functions.GasFunctions.calc_f(Tstart,dgstart,Xgstart);
    parameter Temperature Tstart = T_default;
    parameter SI.Pressure pstart = p_default;
    parameter Real vfracstart = 0.01;

    parameter Real Tinitial = T_pX(pstart,Xfullstart);
  equation
    assert(T >= 273.15 and T <= 573.15, "Temperature T (=" + String(T) + " K = 200 K) is not in the allowed range 273.15 K <= T <= 573.15 K required from medium model \"" + mediumName + "\".");
    Yfull = massToMoleFractions(Xfull,MMX);
    Xg[1:nG-1] = Xfull[1:nG-1]/sum(Xfull[1:nG]);
    sum(Xg) = 1;
    Xl[1:nL-1] = Xfull[1+nG:nG+nL-1]/sum(Xfull[1+nG:nG+nL]);
    vfrac = sum(Xfull[1:nG]);
    sum(Xl) = 1;
    Yg = massToMoleFractions(Xg,MMX[1:nG]);
    Yl = massToMoleFractions(Xl,MMX[1+nG:nG+nL]);
    MM = molarMass(state);
    h =Functions.calc_h(T,dg,Xfull);
    R =Functions.calc_R(Xfull);
    u =Functions.calc_u(T,dg,Xfull);
    d =Functions.calc_d(T,dg,Xfull);
    // connect state with BaseProperties
    state.T = T;
    state.p = p;
    state.X = X;
    state.Xfull = Xfull;
    state.dg = dg;

    p =Functions.GasFunctions.calc_p(T,dg,Xg);
    Xfull = mfull/sum(mfull);
    X = transpose(lambda_mass)*mfull;
    gi =Functions.calc_gi(T, p);
    gim = gi.*MMX;
    m[1:nL-1] =Functions.LiquidFunctions.calc_mfromY(Yl);
    m[nL] = Yl[nL]/IF97.MH2O;
    gamma = Functions.LiquidFunctions.calc_gamma(T,p,Xl);
    a = gamma.*m;
    I =Functions.LiquidFunctions.calc_I(Xl);
    phi = Functions.GasFunctions.calc_phi(T,dg,Xg);
    f = Functions.GasFunctions.calc_f(T,dg,Xg);
    f = exp(logf);
    logreacBase = cat(1,logf,loga);
    gr = nu*gim;
    logK = -gr / Modelica.Constants.R / T;
    K = exp(logK);
    nu*logreacBase = logK;
    a = exp(loga);
    pH = -loga[Hindex-nG]*log10(exp(1));
  end BaseProperties;
//   redeclare function setState_pTX
//     "Return thermodynamic state as function of p, T and composition X"
//     extends Modelica.Icons.Function;
//     input AbsolutePressure p "Pressure";
//     input Temperature T "Temperature";
//     input MassFraction X[:] "Mass fractions";
//     output ThermodynamicState state;
//   protected 
//     MassFraction[nX] Xred = if size(X,1)==nX then X else cat(1,X,{1-sum(X)});
//     MassFraction[nF] Xfull;
//     MassFraction[nG] Xg;
//     SI.Density dg "Gas density";
//   algorithm 
//     Xfull :=X_vfracTXred(p,T,Xred);
//     Xg :=Xfull[1:nG]/sum(Xfull[1:nG]);
//     dg :=dg_TpX(
//       T,
//       p,
//       Xg);
//     state := ThermodynamicState(p=p,T=T,X=Xred,Xfull=Xfull,dg=dg);//,Yred=Yred);
//     annotation (Inline=true, smoothOrder=2);
//   end setState_pTX;

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
    MassFraction[nG] Xg;
    SI.Density dg "Gas density";
  algorithm
    Xfull :=X_phXred(p,h,Xred);
    T :=T_pX(p, Xfull);
    Xg :=Xfull[1:nG]/sum(Xfull[1:nG]);
    dg :=dg_TpX(
      T,
      p,
      Xg);
    state := ThermodynamicState(p=p,T=T,X=Xred,Xfull=Xfull,dg=dg);
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
    MassFraction[nG] Xg;
    SI.Density dg "Gas density";
  algorithm
    Xfull :=X_psXred(p,s,Xred);
    T :=T_pX(p, Xfull);
    Xg :=Xfull[1:nG]/sum(Xfull[1:nG]);
    dg :=dg_TpX(
      T,
      p,
      Xg);
    state := ThermodynamicState(p=p,T=T,X=Xred,Xfull = Xfull,dg=dg);
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
    MassFraction[nG] Xg;
    SI.Density dg "Gas density";
  algorithm
    Xfull :=X_dTXred(d,T,Xred);
    p :=p_TX(T, Xfull);
    Xg :=Xfull[1:nG]/sum(Xfull[1:nG]);
    dg :=dg_TpX(
      T,
      p,
      Xg);
    state := ThermodynamicState(p=p,T=T,X=Xred,Xfull = Xfull,dg = dg);
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
      state.dg,
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
      state.dg,
      state.Xfull);
    annotation(Inline=true,smoothOrder=2);
  end specificEnthalpy;

  redeclare function extends specificInternalEnergy "Return specific internal energy"
    extends Modelica.Icons.Function;
  algorithm
    u :=Functions.calc_u(
      state.T,
      state.dg,
      state.Xfull);
    annotation(Inline=true,smoothOrder=2);
  end specificInternalEnergy;

  redeclare function extends specificEntropy "Return specific entropy"
    extends Modelica.Icons.Function;
  algorithm
    s :=Functions.calc_s(
      state.T,
      state.dg,
      state.Xfull);
    annotation(Inline=true,smoothOrder=2);
  end specificEntropy;

  redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
    extends Modelica.Icons.Function;
  algorithm
    g :=Functions.calc_g(
      state.T,
      state.dg,
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
      state.dg,
      state.Xfull);
    annotation(Inline=true,smoothOrder=1);
  end specificHeatCapacityCp;

  redeclare function extends specificHeatCapacityCv
    "Compute specific heat capacity at constant volume from temperature and gas data"
  algorithm
    cv :=Functions.calc_cv(
      state.T,
      state.dg,
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
end SaturatedLiquid;
