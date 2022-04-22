within ElectrolyteMedia.Media.GasLiquidPhase.Common;
package MixtureLiquid
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

constant Media.GasLiquidPhase.Common.GasInteractionDataRecord
  interactiong=userInterface.interactiong;
  constant LiquidInteractionDataRecord interactionl=userInterface.interactionl;

  constant Integer nL=userInterface.nL;
  constant Integer nLi = nL-1;
  constant Integer nG=userInterface.nG;
  constant Integer nR=userInterface.nR;
  constant Integer nF = nX + nR;

  constant Media.Common.Types.LiquidModel LiquidModel=userInterface.LiquidModel;
  constant Media.Common.Types.GasModel GasModel=userInterface.GasModel;
  constant Real[nR,nF] nu = userInterface.nu;
  constant Real[nR,nF] nu_mass = Functions.Reaction.calc_nu_mass(nu);


  constant UserInterface userInterface annotation (Placement(transformation(extent={{-10,-8},{10,12}})));

constant Real[:,:] lambda=
    Media.Common.Reaction.calc_lambda(
    nu,
    nR,
    nF);
    constant Real[:,:] lambda_mass={{lambda[i,j]/MMX[i] for j in 1:nX} for i in 1:nF}*1000 "Null space of mass based stoichiometry matrix";

  constant MolarMass MH2O = IF97.MH2O "Molar mass of solvent";
  constant SpecificHeatCapacity RH2O = IF97.RH2O;

constant MolarMass[nG + nL] MMX= cat(1,Functions.GasFunctions.calc_MMX(),Functions.LiquidFunctions.calc_MMX());
                                                                                                               // Functions.calc_MMX();

  constant Integer Hindex = max(2,Modelica.Math.BooleanVectors.firstTrueIndex({abs(datal[i].MM-0.001008) == 0 for i in 1:nL-1})+nG);
  constant Integer OHindex = Modelica.Math.BooleanVectors.firstTrueIndex({abs(datal[i].MM-0.017008) == 0 for i in 1:nL-1})+nG;

  redeclare function extends gasConstant "Return gasConstant"
  algorithm
    R :=Functions.calc_R(state.Xfull);
    annotation(Inline = true, smoothOrder = 3);
  end gasConstant;

  redeclare replaceable model extends BaseProperties(
    T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    Xi(each stateSelect=if preferredMediumStates then StateSelect.always else StateSelect.default),
    final standardOrderComponents=true)
    "Base properties (p, d, T, h, u, R, MM, X, and Xfull of aqueous liquid solution"
    MoleFraction[nG+nL] Yfull(start=Yfullstart) "Full mole fraction vector";
    MassFraction[nG+nL] Xfull(start=Xfullstart,each min = 0) "Full mass fraction vector";
    SI.Mass[nG+nL] mfull(start=mfullstart*scale)  "Full amount of substance vector (substitute)";
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
    Real vfrac "Mass based vapor fraction";
    parameter MoleFraction[nG + nL] Yfullstart = Functions.calc_Y_M(Xfullstart, MMX);
  parameter MassFraction[nG + nL] Xfullstart=X_pTXred(
        pstart,
        Tstart,
        Xredstart);
    parameter MassFraction[nX] Xredstart = reference_X;
    parameter Real[nX] Xredstart_ = transpose(lambda_mass)*Xfullstart "for scaling of species moles consistent with maxx invariant vector X";
    parameter Real[nG+nL] mfullstart = Xfullstart/sum(Xredstart_);
    parameter Real[nL-1] mstart_ = Functions.LiquidFunctions.calc_mfromX(Xlstart);
    parameter Real[nL] mstart = cat(1,mstart_,{Ylstart[nL]/MH2O});
    parameter Real[nL] gammastart = Functions.LiquidFunctions.calc_gamma(Tstart,pstart,Xlstart);
    parameter Real[nL] astart = gammastart.*mstart;
    parameter Real[nL] logastart = log(astart);
    final parameter MassFraction[nG] Xgstart = Xfullstart[1:nG]/sum(Xfullstart[1:nG]);
    final parameter MassFraction[nL] Xlstart = Xfullstart[1+nG:nG+nL]/sum(Xfullstart[1+nG:nG+nL]);
    final parameter MoleFraction[nL] Ylstart = Functions.LiquidFunctions.calc_Y(Xlstart);
    final parameter Real[nG] fstart = Functions.GasFunctions.calc_f(Tstart,dgstart,Xgstart);
    parameter SI.Density dgstart = dg_TpX(Tstart,pstart,Xgstart) "Start value for gas density";
    parameter Temperature Tstart = T_default;
    parameter SI.Pressure pstart = p_default;
    parameter Real[:,:] lambda_mass_scaled = lambda_mass/scale;
    parameter Real scale = 100000;//1000000000;//
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
    //     X = transpose(lambda_mass)*mfull;
    X = transpose(lambda_mass_scaled)*mfull;
    // gi =Functions.calc_gi(T, p);
    gi[1:nG] = Functions.GasFunctions.calc_g_i(T);
    gi[1+nG:nG+nL] = Functions.LiquidFunctions.calc_gi(T,p);
    for i in 1:nF loop
      gim[i] = gi[i]*MMX[i];
    end for;
    //     gim = gi.*MMX;
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

  redeclare function setState_pTX
    "Return thermodynamic state as function of p, T and composition X"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:] "Mass fractions";
    output ThermodynamicState state;
protected
    MassFraction[nX] Xred = if size(X,1)==nX then X else cat(1,X,{1-sum(X)});
    MassFraction[nF] Xfull;
    MassFraction[nG] Xg;
    SI.Density dg "Gas density";
  algorithm
    Xfull :=X_pTXred(
      p,
      T,
      Xred);
    Xg :=Xfull[1:nG]/sum(Xfull[1:nG]);
    dg :=dg_TpX(T,p,Xg);
    state := ThermodynamicState(p=p,T=T,X=Xred,Xfull=Xfull,dg=dg);//,Yred=Yred);
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
    MassFraction[nG] Xg;
    SI.Density dg "Gas density";
  algorithm
       (Xfull,T) :=XT_phXred(
         p,
         h,
         Xred);
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
    (Xfull,T) :=XT_psXred(p,s,Xred);
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
    (Xfull,p) :=Xp_dTXred(d,T,Xred);
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
//   redeclare function extends velocityOfSound "Return velocity of sound"
//     extends Modelica.Icons.Function;
//   protected 
//     MolarVolume v = molarVolume(state);
//     MolarMass MM = molarMass(state);
//     Real kappa(unit="1/Pa") = isothermalCompressibility(state);
//   algorithm 
//     a := sqrt(max(0,v/(kappa*MM)));
//     annotation(Inline=true,smoothOrder=2);
//   end velocityOfSound;
// 
//   redeclare function extends isobaricExpansionCoefficient
//     "Returns overall the isobaric expansion coefficient beta"
//   protected 
//     MoleFraction[:] Y = massToMoleFractions(state.X,MMX);
//   algorithm 
//     beta := calc_beta(
//         state.T,
//         state.p,
//         state.X,
//         data,
//         nX);
//     annotation(Inline=true,smoothOrder=2);
//   end isobaricExpansionCoefficient;
// 
//   redeclare function extends isothermalCompressibility
//     "Returns overall the isothermal compressibility factor"
//   protected 
//     MoleFraction[:] Y = massToMoleFractions(state.X,MMX);
//   algorithm 
//     kappa := Functions.calc_kappa(
//         state.T,
//         state.p,
//         state.X,
//         data,
//         nX);
//     annotation(Inline=true,smoothOrder=2);
//   end isothermalCompressibility;
// 
//   redeclare function extends density_derp_T
//     "Returns the partial derivative of density with respect to pressure at constant temperature"
//   algorithm 
//     ddpT := Functions.calc_ddpT(
//         state.T,
//         state.p,
//         state.X,
//         data,
//         nX);
//     annotation(Inline=true,smoothOrder=2);
//   end density_derp_T;
// 
//   redeclare function extends density_derT_p
//     "Returns the partial derivative of density with respect to temperature at constant pressure"
//   algorithm 
//     ddTp := Functions.calc_ddTp(
//         state.T,
//         state.p,
//         state.X,
//         data,
//         nX);
//     annotation(Inline=true,smoothOrder=2);
//   end density_derT_p;
// 
//   function gasMixtureViscosity
//     "Return viscosities of gas mixtures at low pressures (Wilke method)"
//     extends Modelica.Icons.Function;
//     input MoleFraction[:] yi "Mole fractions";
//     input MolarMass[size(yi, 1)] M "Mole masses";
//     input DynamicViscosity[size(yi, 1)] eta "Pure component viscosities";
//     output DynamicViscosity etam "Viscosity of the mixture";
//   protected 
//     Real fi[size(yi, 1),size(yi, 1)];
//   algorithm 
//     for 
//    i in 1:size(eta,1) loop
//    assert(fluidConstants[i].hasDipoleMoment, "Dipole moment for " +
//      fluidConstants[i].chemicalFormula +
//     " not known. Can not compute viscosity.");
//    assert(fluidConstants[i].hasCriticalData, "Critical data for " +
//      fluidConstants[i].chemicalFormula +
//     " not known. Can not compute viscosity.");
//    for 
//      j in 1:size(eta,1) loop
//    if i==1 then
//      fi[i,j] := (1 + (eta[i]/eta[j])^(1/2)*(M[j]/M[i])^(1/4))^2/(8*(1 + M[i]/M[j]))^(1/2);
//    elseif j<i then
//        fi[i,j] := eta[i]/eta[j]*M[j]/M[i]*fi[j,i];
//      else
//        fi[i,j] := (1 + (eta[i]/eta[j])^(1/2)*(M[j]/M[i])^(1/4))^2/(8*(1 + M[i]/M[j]))^(1/2);
//    end if;
//    end for;
//     end for;
//     etam := sum(yi[i]*eta[i]/sum(yi[j]*fi[i, j] for j in 1:size(eta, 1))
//    for i in 1:size(eta, 1));
// 
//     annotation (smoothOrder=2,
//           Documentation(info="<html>
//   
//   <p>
//   Simplification of the kinetic theory (Chapman and Enskog theory)
//   approach neglecting the second-order effects.<br>
//   <br>
//   This equation has been extensively tested (Amdur and Mason, 1958;
//   Bromley and Wilke, 1951; Cheung, 1958; Dahler, 1959; Gandhi and Saxena,
//   1964; Ranz and Brodowsky, 1962; Saxena and Gambhir, 1963a; Strunk, et
//   al., 1964; Vanderslice, et al. 1962; Wright and Gray, 1962). In most
//   cases, only nonpolar mixtures were compared, and very good results
//   obtained. For some systems containing hydrogen as one component, less
//   satisfactory agreement was noted. Wilke's method predicted mixture
//   viscosities that were larger than experimental for the H2-N2 system,
//   but for H2-NH3, it underestimated the viscosities. <br>
//   Gururaja, et al. (1967) found that this method also overpredicted in
//   the H2-O2 case but was quite accurate for the H2-CO2 system. <br>
//   Wilke's approximation has proved reliable even for polar-polar gas
//   mixtures of aliphatic alcohols (Reid and Belenyessy, 1960). The
//   principal reservation appears to lie in those cases where Mi&gt;&gt;Mj
//   and etai&gt;&gt;etaj.<br>
//   </p>
//   
//   </html>"));
//   end gasMixtureViscosity;
// 
  redeclare replaceable function extends dynamicViscosity
    "Return mixture dynamic viscosity"
  algorithm
  eta := 1;//Common.Functions.IF97_R1_Tp.calc_eta(state.p,state.T);
  annotation (smoothOrder=2);
  end dynamicViscosity;
// 
//   function XT_phXred
//     "Liquid dissociation equilibrium calculated via Newton method with homotopy continuation method of step-by step adding solutes to initially pure aqueous solution"
//     input Modelica.SIunits.Pressure p;
//     input SI.SpecificEnthalpy h;
//     input MassFraction[nX] Xred;
//     output Real[nF] X;
//     output Real T;
// 
//   protected 
//     parameter Real eps = 1e-5;
// 
//     SI.Temperature Tmax = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p);
//     SI.Temperature Tmin = 273.15;
// 
//     MoleFraction[nX] Yred;
// 
//     Real[1] f_T;
//     Real[1,1] J_T;
//     SI.Temperature T1;
//     SI.Temperature T2;
//     SI.Temperature[1] Tt_temp;
// 
//     Real[nF] f_X;
//     Real[nF,nF] J_X;
//     MoleFraction[nF] Y1;
//     MoleFraction[nF] Y2;
//   algorithm 
//     //calculate start values for inner-outer algorithm
//     T1 :=exp((log(Tmin) + log(Tmax))/2);
//     X :=X_pTXred(p,T1,Xred);
//     Y1 :=Functions.calc_Y(X,data,nF);
//     Yred :=transpose(lambda)*Y1;
// 
//     //inner (equilibrium) - outer (enthalpy) algorithm
//     while abs(T2/T1-1) > eps loop
//       f_T := Initialization.calc_f.Enthalpy(T1,p,h,X);
//       J_T := Initialization.calc_J.Enthalpy(T1,p,X);
//       Tt_temp:=Initialization.NewtonStep({T1},f_T,J_T,1);
//       assert(sum(Tt_temp)>0,"Error in Newton step for enthalpy calculation",AssertionLevel.error);
//       T2 :=T1;
//       T1 :=Tt_temp[1];
// 
//       while Modelica.Math.Vectors.norm(Y2-Y1) > eps loop
//         Y1 :=Functions.calc_Y(X,data,nF);
//         f_X :=Initialization.calc_f.LE_Tp(
//           Y1,
//           Yred,
//           T1,
//           p);
//         J_X :=Initialization.calc_J.LE_Tp(
//           Y1,
//           T1,
//           p);
//         Y2 := Initialization.NewtonStep(Y1,f_X,J_X,nF);
//         assert(sum(Y2)>0,"Error in Newton step for equilibrium calculation",AssertionLevel.error);
//         Y2 :=Y2/sum(Y2);
//        end while;
//         X :=Functions.calc_X(Y2,data,nF);
//   end while;
// 
//   X :=X_pTXred(p,T1,Xred);
//   T :=T1;
//   end XT_phXred;
// 
//   function XT_psXred
//     "Liquid dissociation equilibrium calculated via Newton method with homotopy continuation method of step-by step adding solutes to initially pure aqueous solution"
//     input Modelica.SIunits.Pressure p;
//     input SI.SpecificEntropy s;
//     input MassFraction[nX] Xred;
//     output Real[nF] X;
//     output Real T;
// 
//   protected 
//     parameter Real eps = 1e-5;
// 
//     SI.Temperature Tmax = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p);
//     SI.Temperature Tmin = 273.15;
// 
//     MoleFraction[nX] Yred;
// 
//     Real[1] f_T;
//     Real[1,1] J_T;
//     SI.Temperature T1;
//     SI.Temperature T2;
//     SI.Temperature[1] Tt_temp;
// 
//     Real[nF] f_X;
//     Real[nF,nF] J_X;
//     MoleFraction[nF] Y1;
//     MoleFraction[nF] Y2;
//   algorithm 
//     //calculate start values for inner-outer algorithm
//     T1 :=exp((log(Tmin) + log(Tmax))/2);
//     X :=X_pTXred(p,T1,Xred);
//     Y1 :=Functions.calc_Y(X,data,nF);
//     Yred :=transpose(lambda)*Y1;
// 
//     //inner (equilibrium) - outer (enthalpy) algorithm
//     while abs(T2/T1-1) > eps loop
//       f_T := Initialization.calc_f.Entropy(T1,p,s,X);
//       J_T := Initialization.calc_J.Entropy(T1,p,X);
//       Tt_temp:=Initialization.NewtonStep({T1},f_T,J_T,1);
//       assert(sum(Tt_temp)>0,"Error in Newton step for enthalpy calculation",AssertionLevel.error);
//       T2 :=T1;
//       T1 :=Tt_temp[1];
// 
//       while Modelica.Math.Vectors.norm(Y2-Y1) > eps loop
//         Y1 :=Functions.calc_Y(X,data,nF);
//         f_X :=Initialization.calc_f.LE_Tp(
//           Y1,
//           Yred,
//           T1,
//           p);
//         J_X :=Initialization.calc_J.LE_Tp(
//           Y1,
//           T1,
//           p);
//         Y2 := Initialization.NewtonStep(Y1,f_X,J_X,nF);
//         assert(sum(Y2)>0,"Error in Newton step for equilibrium calculation",AssertionLevel.error);
//         Y2 :=Y2/sum(Y2);
//        end while;
//         X :=Functions.calc_X(Y2,data,nF);
//   end while;
// 
//   X :=X_pTXred(p,T1,Xred);
//   T :=T1;
//   end XT_psXred;
// 
//   function Xp_dTXred
//     "Liquid dissociation equilibrium calculated via Newton method with homotopy continuation method of step-by step adding solutes to initially pure aqueous solution"
//     input SI.Density d;
//     input SI.Temperature T;
//     input MassFraction[nX] Xred;
//     output Real[nF] X;
//     output SI.Pressure p;
// 
//   protected 
//     parameter Real eps = 1e-5;
// 
//     SI.Pressure pmin = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.psat(T);
//     SI.Pressure pmax = 1000e5;
// 
//     MoleFraction[nX] Yred;
// 
//     Real[1] f_p;
//     Real[1,1] J_p;
//     SI.Pressure p1;
//     SI.Pressure p2;
//     SI.Temperature[1] p_temp;
// 
//     Real[nF] f_X;
//     Real[nF,nF] J_X;
//     MoleFraction[nF] Y1;
//     MoleFraction[nF] Y2;
//   algorithm 
//     //calculate start values for inner-outer algorithm
//     p1 :=exp((log(pmin) + log(pmax))/2);
//     X :=X_pTXred(p1,T,Xred);
//     Y1 :=Functions.calc_Y(X,data,nF);
//     Yred :=transpose(lambda)*Y1;
// 
//     //inner (equilibrium) - outer (enthalpy) algorithm
//     while abs(p2/p1-1) > eps loop
//       f_p := Initialization.calc_f.Density(T,p1,X,d);
//       J_p := Initialization.calc_J.Density(T,p1,X);
//       p_temp:=Initialization.NewtonStep({p1},f_p,J_p,1);
//       assert(sum(p_temp)>0,"Error in Newton step for enthalpy calculation",AssertionLevel.error);
//       p2 :=p1;
//       p1 :=p_temp[1];
// 
//       while Modelica.Math.Vectors.norm(Y2-Y1) > eps loop
//         Y1 :=Functions.calc_Y(X,data,nF);
//         f_X :=Initialization.calc_f.LE_Tp(
//           Y1,
//           Yred,
//           T,
//           p1);
//         J_X :=Initialization.calc_J.LE_Tp(
//           Y1,
//           T,
//           p1);
//         Y2 := Initialization.NewtonStep(Y1,f_X,J_X,nF);
//         assert(sum(Y2)>0,"Error in Newton step for equilibrium calculation",AssertionLevel.error);
//         Y2 :=Y2/sum(Y2);
//        end while;
//         X :=Functions.calc_X(Y2,data,nF);
//   end while;
// 
//   X :=X_pTXred(p1,T,Xred);
//   p :=p1;
//   end Xp_dTXred;

end MixtureLiquid;
