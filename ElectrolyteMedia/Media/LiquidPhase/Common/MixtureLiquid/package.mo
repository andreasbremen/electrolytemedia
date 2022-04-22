within ElectrolyteMedia.Media.LiquidPhase.Common;
package MixtureLiquid "Medium model of a mixture of liquids based on Modelica library"
  import Modelica.Math;
  import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;

  extends Modelica.Media.Interfaces.PartialMixtureMedium(
     ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX,
     fixedX=false,
     reducedX = true,
     singleState=false,
     nX = nS-nR,
     nXi = nX-1,
     substanceNames = userInterface.substanceNames,
     mediumName = userInterface.mediumName,
     T_default = userInterface.Tstart,
     p_default = userInterface.pstart,
     reference_X = userInterface.refX,
     Density(start=10, nominal=10),
     AbsolutePressure(start=10e5, nominal=10e5),
     Temperature(min=273.15, max=573.15, start=300, nominal=300));
     //substanceNames=data[:].name,
     // fixedX=true,
     // reducedX = false,
     // singleState=false,
//      nXi = nS-nR-1,
//      
//      SpecificEnthalpy(start=if Functions.referenceChoice==ReferenceEnthalpy.ZeroAt25C then 3e5 else 
//         if Functions.referenceChoice==ReferenceEnthalpy.UserDefined then Functions.h_offset else 0, nominal=1.0e5),

  redeclare record extends ThermodynamicState "Thermodynamic state variables"
    Real[nF] Xfull;
  end ThermodynamicState;
  constant UserInterface userInterface "User interface to pass data"  annotation (Placement(transformation(extent={{-10,-8},{10,12}})));
  constant DataRecord[nLi] data = userInterface.data "Solutes data record from user interface";
  constant Media.Common.Types.LiquidModel LiquidModel = userInterface.LiquidModel "Liquid gibbs excess model from user interface";
  constant LiquidInteractionDataRecord interaction=userInterface.interaction "Bromley binary interaction data record from user interface";

  constant Integer nR=userInterface.nR "Number of gas-liquid and dissociation equilibria from user interface";
  constant Integer nF = nX + nR "Number of all species in all phases";
  constant Integer nL = userInterface.nL "Number of liquid species (solutes+solvent) from user interface";
  constant Integer nLi = nL-1 "Number of liquid species (solutes+solvent) -1";
  constant SI.MassFraction[nL] reference_Xfull= userInterface.refXfull "Reference mass fraction vector of all species in the system";
                              //max(1,nF-1);
  constant Real[nR,nF] nu = userInterface.nu "Stoichiometry matrix of gas-liquid and dissociation equilibria from user interface";
  constant Real[nR,nF] nu_mass = Functions.Reaction.calc_nu_mass(nu) "Mass based stoichiometry matrix of gas-liquid and dissociation equilibria from user interface";

  constant Real[:,:] lambda=
      Media.Common.Reaction.calc_lambda(
      nu,
      nR,
      nF)
         "Null space of stoichiometry matrix";
  constant Real[:,:] lambda_mass= {{lambda[i,j]/MMX[i] for j in 1:nX} for i in 1:nL}*1000 "Null space of mass based stoichiometry matrix";

  constant MolarMass MH2O = IF97.MH2O "Molar mass of solvent";
  constant SpecificHeatCapacity RH2O = IF97.RH2O "Gas constant of H2O";

  constant MolarMass[nL] MMX= Functions.calc_MMX() "Vector containing molar masses of all species";

  constant Integer Hindex(min=1) = max(1,Modelica.Math.BooleanVectors.firstTrueIndex({abs(data[i].MM-0.001008) < Modelica.Constants.eps for i in 1:nF-1}))        "Index of H+ ion";
  constant Integer OHindex(min=1) = Modelica.Math.BooleanVectors.firstTrueIndex({abs(data[i].MM-0.017008) < Modelica.Constants.eps for i in 1:nF-1})         "Index of OH- ion";

  redeclare function extends gasConstant "Return gasConstant"
  algorithm
    R := sum(data.R.*state.Xfull[1:nF-1]) + RH2O*state.Xfull[nF];
    annotation(Inline = true, smoothOrder = 3);
  end gasConstant;

  redeclare replaceable model extends BaseProperties(
    T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
    Xi(each stateSelect=if preferredMediumStates then StateSelect.always else StateSelect.default),
    final standardOrderComponents=true)
    "Base properties (p, d, T, h, u, R, MM, X, and Xfull of aqueous liquid solution"

    MoleFraction[nF] Yfull(start=Yfullstart) "Full mole fraction vector";
    MassFraction[nF] Xfull(start=Xfullstart) "Full mass fraction vector";
    SI.Mass[nF] mfull(start=mfullstart)  "Full amount of substance vector (substitute)";
    SI.MolarEnergy[nF] gi;
    SI.MolarEnergy[nR] gr "Gibbs free energy of reaction at T and p ";
    Real[nR] K "Equilibrium constant";
    Real[nR] logK "Logarithmic equilibrium constant";
    Real[nF] m "Molality";
    Real[nF] a(start=astart) "Activity";
    Real[nF] gamma "Activity coefficient";
    Real pH(start = -logastart[Hindex]*log10(exp(1)));
    Real[nF] loga(start=logastart);
    Real I "Molality based ionic strength";
    Real[nX] X_ "Substitute for X to ensure non-reactive species to be > 0";
    parameter MoleFraction[nF] Yfullstart= Functions.calc_Y(Xfullstart);
    parameter MassFraction[nF] Xfullstart = X_pTXred(pstart,Tstart,Xredstart);
    parameter MassFraction[nX] Xredstart = reference_X;
    parameter Real[nX] Xredstart_ = transpose(lambda_mass)*Xfullstart "for scaling of species moles consistent with maxx invariant vector X";
    parameter Real[nF] mstart = cat(1,Functions.calc_mfromX(Xfullstart),{Yfullstart[end]/MH2O});
    parameter Real[nF] gammastart = Functions.calc_gamma(Tstart,pstart,Xfullstart);
    parameter Real[nF] astart = gammastart.*mstart;
    parameter Real[nF] logastart = log(astart);
    parameter Real[nF] mfullstart = Xfullstart/sum(Xredstart_);
    parameter Temperature Tstart = T_default;
    parameter SI.Pressure pstart = p_default;

  equation

    assert(T >= 273.15 and T <= 573.15, "Temperature T (="
                  + String(T) + " K = 200 K) is not in the allowed range
                273.15 K <= T <= 573.15 K required from medium model \""   + mediumName + "\".");

    Yfull = massToMoleFractions(Xfull,MMX);
    MM = molarMass(state);
    h =Functions.calc_h(T,p,Xfull);
    R = cat(1,data[:].R,{RH2O})*Xfull;
    u =Functions.calc_u(T,p,Xfull);
    d =Functions.calc_d(T,p,Xfull);
    // connect state with BaseProperties
    state.T = T;
    state.p = p;
    state.X = X;
    state.Xfull = Xfull;

    Xfull = mfull/sum(mfull);
    X_ = transpose(lambda_mass)*mfull;
    X_ = {max(1e-25,X[i]) for i in 1:nX};

    gi =Functions.calc_gi(T, p) .* MMX;
    m[1:nF-1] =Functions.calc_mfromY(Yfull);
    m[nF] = Yfull[nF]/MH2O;
    gamma =Functions.calc_gamma(T,p,Xfull);
    a = gamma.*m;
    gr = nu*gi;
    logK = -gr / Modelica.Constants.R / T;
    K = exp(logK);
    nu*loga = logK;
    a = exp(loga);
    pH = -loga[Hindex]*log10(exp(1));
    I = Functions.calc_I(Xfull);
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
  algorithm
    Xfull :=X_pTXred(p,T,Xred);
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
    (Xfull,p) :=Xp_dTXred(
      d,
      T,
      Xred);
    state := ThermodynamicState(p=p,T=T,X=Xred,Xfull = Xfull);
    annotation (Inline=true, smoothOrder=2);
  end setState_dTX;

  redeclare function extends pressure "Return pressure"
  algorithm
    p := state.p;
    annotation (Inline=true, smoothOrder=2);
  end pressure;

  redeclare function extends temperature "Return temperature"
  algorithm
    T := state.T;
    annotation(Inline=true,smoothOrder=2);
  end temperature;

  redeclare function extends density "Return density"
  algorithm
    d :=Functions.calc_d(
        state.T,
        state.p,
        state.Xfull);
    annotation(Inline = true, smoothOrder = 3);
  end density;

  redeclare function extends molarMass "Return molar mass of mixture"
  algorithm
    MM := 1/(sum(state.Xfull[j]/data[j].MM for j in 1:nF-1) + state.Xfull[nF]/MH2O);
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
//    redeclare function extends isobaricExpansionCoefficient
//      "Returns overall the isobaric expansion coefficient beta"
//    algorithm 
//      beta := Functions.calc_beta(
//          state.T,
//          state.p,
//          state.Xfull);
//      annotation(Inline=true,smoothOrder=2);
//    end isobaricExpansionCoefficient;
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
  eta :=
      Functions.IF97_R1_Tp.calc_eta(
      state.p, state.T);
  annotation (smoothOrder=2);
  end dynamicViscosity;
end MixtureLiquid;
