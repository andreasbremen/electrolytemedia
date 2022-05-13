within ElectrolyteMedia.Media.SolidLiquidPhase;
package Common

  constant Modelica.SIunits.Pressure pref = 1e5;
  constant Modelica.SIunits.Temperature Tref = 298.15;
  constant Real eta = 1.66027e5 "constant for w calculations";
  constant Modelica.SIunits.Pressure Psi = 2600e5 "solvent characerstic constant in Pa";
  constant Modelica.SIunits.Temperature Theta = 228 "solvent characteristic constant in K";
  constant Real Zref = -1.278034682000000E-002 "Z at reference state";
  constant Real Yref = -5.798650444000000E-005 "Y at reference state";

  record DataRecordL
    extends LiquidPhase.Common.DataRecord;
  end DataRecordL;

  record DataRecordS
    extends SolidPhase.Common.DataRecord;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end DataRecordS;

  record LiquidInteractionDataRecord
    extends LiquidPhase.Common.LiquidInteractionDataRecord;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end LiquidInteractionDataRecord;

  record ConstDataRecord
    constant Real eta = 1.66027e5 "constant for w calculations";
    constant Modelica.SIunits.Pressure Psi = 2600e5 "solvent characerstic constant in Pa";
    constant Modelica.SIunits.Temperature Theta = 228 "solvent characteristic constant in K";
    constant Real Z_ref = -1.278034682000000E-002 "Z at reference state";
    constant Real Y_ref = -5.798650444000000E-005 "Y at reference state";
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end ConstDataRecord;

  record IF97 "record for IF97 parameters of region 1; pT explicit"
    extends Modelica.Media.Water.IF97_Utilities.BaseIF97.data(MH2O = 0.018016);
    constant SI.Temperature Ttriple=273.16 "The triple point temperature";
    constant SI.Pressure ptriple=611.657 "The triple point pressure";
    constant SI.Density dltriple=999.792520031617642
      "The triple point liquid density";
    constant SI.Density dvtriple=0.485457572477861372e-2
      "The triple point vapour density";
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end IF97;

  record StdRefH2O
    constant SI.SpecificEnthalpy h_tr = -287720.759/0.0180158 "data from Helgeson1974";
    constant SI.SpecificEnergy g_tr = -235512.710/0.0180158 "data from Helgeson1974";
    constant SI.SpecificEnergy u_tr = -284027.643/0.0180158 "data from Helgeson1974";
    constant SI.SpecificEnergy a_tr = -231855.624/0.0180158 "data from Helgeson1974";
    constant SI.SpecificEntropy s_tr = 63.31262/0.0180158 "data from Helgeson1974";
    constant SI.Temperature T_tr = 273.16 "triple point temperature";
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end StdRefH2O;

  record Born
    constant Real a_i[:] = {14.70333593, 212.8462733, -115.4445173, 19.55210915, -83.3034798, 32.13240048, -6.694098645, -37.86202045, 68.87359646, -27.29401652};
    constant SI.Temperature T_eps = 298.15;
    constant SI.Density rho_eps = 1000;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Born;

  record Solvent
    constant Real a[:] = {-2.037662, 5.747e-3, -6.557892e-6};
    constant Real b[:] = {6.107361, -1.074377e-2, 1.268348e-5};
    constant Modelica.SIunits.Density rho_0 = 1000;
    constant Real a_f[:] = {3.66666e-16, -1.504956e-10, 5.01799e-14};
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Solvent;

  record GibbsDerivs
    extends Modelica.Media.Common.GibbsDerivs;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end GibbsDerivs;

  record ThermoProperties_pT
    extends Modelica.Media.Common.ThermoFluidSpecial.ThermoProperties_pT(d(min=0));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end ThermoProperties_pT;

  package SolutesData
    extends Media.LiquidPhase.Common.SolutesData;

  end SolutesData;

  package MixtureSolutesData
    extends Media.LiquidPhase.Common.MixtureSolutesData;

  end MixtureSolutesData;

  package SingleSolidsData
    extends SolidPhase.Common.SolidData;

  end SingleSolidsData;

  package Functions

    constant Integer nSfunfun(min=1)=1;
    constant DataRecordS[nSfunfun] dataSfunfun;

    constant Integer nLfunfun=2;
    final constant Integer nLifunfun = nLfunfun-1;
    constant Common.DataRecordL[nLfunfun - 1] dataLfunfun;
    constant LiquidInteractionDataRecord interactionLfunfun;
    constant Media.Common.Types.LiquidModel LiquidModelfunfun;

    package Reaction "Package containing all functions regarding dissociation reactions"

      function nullSpace "Return the orthonormal nullspace of a matrix"
        extends Modelica.Icons.Function;

        input Real[:,:] nu;
        output Real[size(nu,2),size(nu,2)-size(nu,1)] nullSpace;
      protected
        Integer nullity;
      algorithm

        (nullSpace,nullity) :=Modelica.Math.Matrices.nullSpace(nu);

      //   input Real A[nL, nY] "Input matrix";
      //   input Integer nL;
      //   input Integer nY;
      //   output Real Z[size(A, 2), nY-nL] "Orthonormal nullspace of matrix A";
      //   //output Integer nullity "Nullity, i.e., the dimension of the nullspace";
      //
      // protected
      //   Real V[size(A, 2), size(A, 2)] "Right orthogonal matrix";
      //   Real sigma[min(size(A, 1), size(A, 2))] "singular values";
      //   Integer rank "rank of matrix A";
      //   Real eps "tolerance for rank determination";
      //   Integer n=min(size(A, 1), size(A, 2));
      //   Integer i=n;
      //   Integer nullity;
      //
      // algorithm
      //   (sigma,,V) := Modelica.Math.Matrices.singularValues(A);
      //   V := transpose(V);
      //   // rank computation
      //   eps := max(size(A, 1), size(A, 2))*max(sigma)*Modelica.Constants.eps;
      //   rank := 0;
      //   if n > 0 then
      //     while i > 0 loop
      //       if sigma[i] > eps then
      //         rank := i;
      //         i := 0;
      //       end if;
      //       i := i - 1;
      //     end while;
      //   end if;
      //   Z := V[:, rank + 1:size(A, 2)];
      //   // nullspace computation
      //   nullity := size(A, 2) - rank;
      //   // nullity
      end nullSpace;

      function calc_nu_mass
        input Real[:,:] nu;

        output Real[size(nu,1),size(nu,2)] nu_mass;
      protected
        SI.MolarMass[size(nu,2)] MMX = calc_MMX();
      algorithm
        for i in 1:size(nu,1) loop
          nu_mass[i,:] :=nu[i, :] .* MMX[:];
        end for;

      //   if GasModelfunfun == Media.Common.Types.GasModel.Ideal then
      //     for i in 1:nSfunfun loop
      //       nu_mass[:,i] :=nu[:, i] * dataIGfunfun[i].MM;
      //     end for;
      //   elseif GasModelfunfun == Media.Common.Types.GasModel.PengRobinson then
      //     for i in 1:nSfunfun loop
      //       nu_mass[:,i] :=nu[:,i] * dataPRfunfun[i].MM;
      //     end for;
      //   end if;
      //   for i in 1:nSfunfun:nSfunfun+nLfunfun-1 loop
      //     nu_mass[:,i] :=nu[:,i] * dataLfunfun[i].MM;
      //   end for;
      //   nu_mass[:,nSfunfun+nLfunfun] :=nu[:, nSfunfun + nLfunfun] * IF97.MH2O;
      end calc_nu_mass;
    end Reaction;

    function calc_h
      "calculates specific enthalpy of gas, solid mixture at T and p"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output SI.SpecificEnthalpy h;

    protected
      SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
      SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
      SI.SpecificEnthalpy hs = SolidFunctions.calc_h(T,p,Xs);
      SI.SpecificEnthalpy hl = LiquidFunctions.calc_h(T,p,Xl);
    algorithm

      h:=sum(X[1:nSfunfun])*hs + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*hl;

      annotation(Inline=true,smoothOrder=5);
    end calc_h;

    function calc_s
      "calculates specific entropy of gas, solid mixture at T and p"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output SI.SpecificEntropy s;

    protected
      SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
      SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
      SI.SpecificEntropy ss = SolidFunctions.calc_s(T,Xs);
      SI.SpecificEntropy sl = LiquidFunctions.calc_s(T,p,Xl);
    algorithm

      s:=sum(X[1:nSfunfun])*ss + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*sl;

      annotation(smoothOrder=5);
    end calc_s;

    function calc_u
      "calculates specific internal energy of gas, solid mixture at T and p"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output SI.SpecificInternalEnergy u;

    protected
      SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
      SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
      SI.SpecificInternalEnergy us = SolidFunctions.calc_u(T,p,Xs);
      SI.SpecificInternalEnergy ul = LiquidFunctions.calc_u(T,p,Xl);
    algorithm

      u:=sum(X[1:nSfunfun])*us + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*ul;

      annotation(Inline=true,smoothOrder=5);
    end calc_u;

    function calc_d
      "calculates density of gas, solid mixture at T and p"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output SI.Density d;
    protected
      SI.SpecificVolume v;
    algorithm

      v :=calc_v(T,p,X);
      d :=1/v;

      annotation(Inline=true,smoothOrder=5);
    end calc_d;

    function calc_cp
      "calculates specific heat capacity cp of solid mixture at T and p"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output SI.SpecificHeatCapacityAtConstantPressure cp;

    protected
      SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
      SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
      SI.SpecificHeatCapacityAtConstantPressure cps = SolidFunctions.calc_cp(T,Xs);
      SI.SpecificHeatCapacityAtConstantPressure cpl = LiquidFunctions.calc_cp(T,p,Xl);
    algorithm

      cp:=sum(X[1:nSfunfun])*cps + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*cpl;

    end calc_cp;

    function calc_cv
      "calculates specific heat capacity cv of solid mixture at T and p"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output SI.SpecificHeatCapacityAtConstantVolume cv;

    protected
      SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
      SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
      SI.SpecificHeatCapacityAtConstantVolume cvs = SolidFunctions.calc_cp(T,Xs);
      SI.SpecificHeatCapacityAtConstantVolume cvl = LiquidFunctions.calc_cv(T,p,Xl);
    algorithm

      cv:=sum(X[1:nSfunfun])*cvs + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*cvl;

    end calc_cv;

    function calc_g
      "calculates specific Gibbs free energy of gas, solid mixture at T and p"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output SI.SpecificGibbsFreeEnergy g;

    protected
      SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
      SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
      SI.SpecificGibbsFreeEnergy gs = SolidFunctions.calc_g(T,p,Xs);
      SI.SpecificGibbsFreeEnergy gl = LiquidFunctions.calc_g(T,p,Xl);
    algorithm

      g:=sum(X[1:nSfunfun])*gs + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*gl;
    annotation(smoothOrder=5);
    end calc_g;

    function calc_gi
      "calculates partial Gibbs free energies of solutes and solvent at T and p"
      input SI.Temperature T;
      input SI.Pressure p;

      output SI.SpecificEnergy[nSfunfun+nLfunfun] gi;
    algorithm
       gi[1:nSfunfun] :=SolidFunctions.calc_g_i(T,p);
       gi[1+nSfunfun:nSfunfun+nLfunfun] :=LiquidFunctions.calc_gi(T,p);

      annotation(smoothOrder=5);
    end calc_gi;

    function calc_ddpT
      "calculates pressure derivative of density of solid mixture at T and p"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output SI.DerDensityByPressure ddpT;

    protected
      SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
      SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
      SI.DerDensityByPressure dsdpT;
      SI.DerDensityByPressure dldpT;
      Real vsdpT;
      Real vldpT;
      Real vdpT;
      SI.Density ds;
      SI.Density dl;
      SI.Density d;
    algorithm
      ds :=SolidFunctions.calc_d(T,p,Xs);
      dl :=LiquidFunctions.calc_d(T,p,Xl);

      dsdpT :=SolidFunctions.calc_der_d_p(
        T,
        p,
        Xs);
      dldpT :=LiquidFunctions.calc_der_d_p(
        T,
        p,
        Xl);

      vsdpT :=-1/ds^2*dsdpT;
      vldpT :=-1/dl^2*dldpT;

      vdpT :=sum(X[1:nSfunfun])*vsdpT + sum(X[1 + nSfunfun:nSfunfun + nLfunfun])*vldpT;

      d :=calc_d(T, p, X);

      ddpT :=-vdpT*d^2;

    end calc_ddpT;

    function calc_v
      "calculates specific volume of solid mixture at T and p"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output SI.SpecificVolume v;

    protected
      SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
      SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
      SI.SpecificVolume vs = SolidFunctions.calc_v(T,p,Xs);
      SI.SpecificVolume vl = LiquidFunctions.calc_v(T,p,Xl);
    algorithm

      v:=sum(X[1:nSfunfun])*vs + sum(X[1+nSfunfun:nSfunfun+nLfunfun])*vl;
      annotation(smoothOrder=5);
    end calc_v;

    function calc_vdpT
      "calculates pressure derivative of specific volume of solid mixture at T and p"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output Real vdpT;

    protected
      SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
      SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
      SI.DerDensityByPressure dsdpT;
      SI.DerDensityByPressure dldpT;
      Real vsdpT;
      Real vldpT;
      SI.Density ds;
      SI.Density dl;
    algorithm
      ds :=SolidFunctions.calc_d(T,p,Xs);
      dl :=LiquidFunctions.calc_d(T,p,Xl);

      dsdpT :=SolidPhase.Common.Functions.calc_der_d_p(
        T,
        p,
        Xs);
      dldpT :=LiquidPhase.Common.Functions.calc_der_d_p(
        T,
        p,
        Xl);

      vsdpT :=-1/ds^2*dsdpT;
      vldpT :=-1/dl^2*dldpT;

      vdpT :=sum(X[1:nSfunfun])*vsdpT + sum(X[1 + nSfunfun:nSfunfun + nLfunfun])*vldpT;

    end calc_vdpT;

    function calc_MMX
      output SI.MolarMass[nSfunfun+nLfunfun] MMX;

    protected
      SI.MolarMass[nSfunfun] MMXs = SolidFunctions.calc_MMX();
      SI.MolarMass[nLfunfun] MMXl = LiquidFunctions.calc_MMX();
    algorithm
      MMX[1:nSfunfun] :=MMXs;// :=cat(1,MMXs,MMXl);
      MMX[1+nSfunfun:nSfunfun+nLfunfun] :=MMXl;

    end calc_MMX;

    function calc_R
      input Modelica.Media.Interfaces.Types.MassFraction[nSfunfun + nLfunfun] X;
      output SI.SpecificEntropy R;
    protected
      SI.SpecificEntropy Rs = SolidFunctions.calc_R(X[1:nSfunfun]);
      SI.SpecificEntropy Rl = LiquidFunctions.calc_R(X[1+nSfunfun:nSfunfun+nLfunfun]);
    algorithm
      R :=Rs + Rl;
      annotation(Inline=true,smoothOrder=5);
    end calc_R;

    function calc_MM
      "calculates molar mass of gas, solid mixture at T and p"
      input SI.MassFraction[nSfunfun+nLfunfun] X;

      output SI.MolarMass MM;

    protected
      SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
      SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
      SI.MolarMass MMs = SolidFunctions.calc_MM(Xs);
      SI.MolarMass MMl = LiquidFunctions.calc_MM(Xl);
    algorithm
        MM :=1/(sum(X[1:nSfunfun])/MMs + sum(X[1+nSfunfun:nSfunfun+nLfunfun])/MMl);

    end calc_MM;

    function calc_Y_M "calculates mole fraction vector of liquid species"
      input SI.MassFraction[:] X;
      input SI.MolarMass[size(X,1)] M;
      output SI.MoleFraction[size(X,1)] Y;
    protected
      Real[size(X,1)] X_(unit="mol/kg");
    algorithm
      for i in 1:size(X,1) loop
        X_[i] :=X[i]/M[i];
      end for;
      Y :=X_/sum(X_);

      annotation(smoothOrder=5);
    end calc_Y_M;

    function calc_X_M "calculates mass fraction vector of liquid species"
      input SI.MassFraction[:] Y;
      input SI.MolarMass[size(Y,1)] M;
      output SI.MoleFraction[size(Y,1)] X;
    protected
      Real[size(Y,1)] Y_(unit="kg/mol");
    algorithm
      for i in 1:size(Y,1) loop
        Y_[i] :=Y[i]*M[i];
      end for;
      X :=Y_/sum(Y_);

      annotation(smoothOrder=5);
    end calc_X_M;

    package LiquidFunctions
      extends LiquidPhase.Common.Functions(nLfun=nLfunfun,  datafun=dataLfunfun,  LiquidModelfun=LiquidModelfunfun, interactionfun = interactionLfunfun);

    end LiquidFunctions;

    package SolidFunctions
      extends SolidPhase.Common.Functions(nSfun = nSfunfun, datafun=dataSfunfun);

    end SolidFunctions;

    function calc_Yfull "calculates total mole fraction vector"
      input SI.MassFraction [nLfunfun+nSfunfun] X;
      output SI.MoleFraction[nLfunfun+nSfunfun] Y;
    protected
      Real[nLfunfun+nSfunfun] X_(unit="mol/kg");
    algorithm
     for i in nSfunfun+1:nSfunfun+nLfunfun-1 loop
        X_[i] :=X[i]/dataLfunfun[i-nSfunfun].MM;
     end for;
     for i in 1:nSfunfun loop
        X_[i] :=X[i]/ dataSfunfun[i].MM;//0.05844;//
     end for;
      X_[nSfunfun+nLfunfun] :=X[nSfunfun+nLfunfun]/IF97.MH2O;
      Y :=X_/sum(X_);
    end calc_Yfull;

    function calc_Xfull
      "calculates full mass fraction vector from full mole fraction vector"
      input SI.MoleFraction[nLfunfun+nSfunfun] Y;
      output SI.MassFraction [nLfunfun+nSfunfun] X;

    protected
      Real[nLfunfun+nSfunfun] Y_(unit="mol/kg");

    algorithm
     for i in nSfunfun+1:nSfunfun+nLfunfun-1 loop
        Y_[i] :=Y[i]*dataLfunfun[i-nSfunfun].MM;
     end for;

     for i in 1:nSfunfun loop
        Y_[i] :=Y[i]* dataSfunfun[i].MM;//0.05844;//
     end for;
      Y_[nSfunfun+nLfunfun] :=Y[nSfunfun+nLfunfun]*IF97.MH2O;
      X :=Y_/sum(Y_);

    end calc_Xfull;

    function calc_Xs "calculates mass fraction vector of solid phase"
      input SI.MassFraction[nSfunfun+nLfunfun] X;
      output SI.MassFraction[nSfunfun] Xs;
    algorithm
      if sum(X[1:nSfunfun]) > 0 then
        Xs :=X[1:nSfunfun]/sum(X[1:nSfunfun]);
      else
        Xs :=fill(1/nSfunfun, nSfunfun);
      end if;

      annotation(smoothOrder=5);
    end calc_Xs;

    function calc_Xl "calculates mass fraction vector of liquid phase"
      input SI.MassFraction[nSfunfun+nLfunfun] X;
      output SI.MassFraction[nLfunfun] Xl;
    algorithm
      if sum(X[1+nSfunfun:nSfunfun+nLfunfun]) > 0 then
        Xl :=X[1+nSfunfun:nSfunfun+nLfunfun]/sum(X[1+nSfunfun:nSfunfun+nLfunfun]);
      else
        Xl :=fill(1/nLfunfun, nLfunfun);
      end if;

      annotation(smoothOrder=5);
    end calc_Xl;

    function calc_Phil "calculates liquid volume fraction"
      input SI.Temperature T;
      input SI.Pressure p;
      input SI.MassFraction[nSfunfun+nLfunfun] X;
      output SI.VolumeFraction Phil;
    protected
      SI.MassFraction lfrac = sum(X[1+nSfunfun:nSfunfun+nLfunfun]);
      SI.MassFraction[nSfunfun] Xs = X[1:nSfunfun]/sum(X[1:nSfunfun]);
      SI.MassFraction[nLfunfun] Xl = X[1+nSfunfun:nSfunfun+nLfunfun]/sum(X[1+nSfunfun:nSfunfun+nLfunfun]);
      SI.SpecificVolume vl = LiquidFunctions.calc_v(T,p,Xl);
      SI.SpecificVolume vs = SolidFunctions.calc_v(T,p,Xs);
    algorithm
      Phil :=lfrac*vl/(lfrac*vl + (1 - lfrac)*vs);

      annotation(Inline=true,smoothOrder=5);
    end calc_Phil;
  end Functions;

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

    record UserInterface
      extends Modelica.Icons.Record;

      //Solid phase
      parameter Integer ns = 1 "Number of solid species" annotation(Dialog(tab="Solid phase"));

      String[:] solidSubstanceNames = {"solid_1"} "Solid species names" annotation(Dialog(tab="Solid phase"));

      DataRecordS[ns] datas "Species data record" annotation(choicesAllMatching=true, Dialog(tab="Solid phase"));

      //Liquid phase
      parameter Integer nL = 2 "Number of liquid species (solutes+solvent)" annotation(Evaluate = true, Dialog(tab="Liquid phase"));

      Media.Common.Types.LiquidModel LiquidModel "Liquid Gibbs excess model" annotation(choicesAllMatching=true, Dialog(tab="Liquid phase"));

      String[:] liquidSubstanceNames = {"solute","solvent"} "Liquid species names (solutes+solvent)" annotation(Evaluate = true,Dialog(tab="Liquid phase"));

      DataRecordL[nL-1] datal "Solutes data record" annotation(choicesAllMatching=true, Dialog(tab="Liquid phase"));

      LiquidInteractionDataRecord interactionL "Binary interaction data record" annotation (choicesAllMatching=true, Dialog(tab="Liquid phase"));

      //Reaction
      parameter Integer nR = 1 "Number of solid-liquid and dissociation equilibria" annotation(Dialog(tab="Reaction"));

      parameter Real[nR,ns+nL] nu = cat(2,identity(nR),-ones(nR,ns+nL-nR)) "Stoichiometry matrix of solid-liquid and dissociation equilibria" annotation(Dialog(tab="Reaction"));

      //General
      String mediumName = "mediumName" "Name of medium" annotation(Dialog(tab="General"));

      SI.Temperature Tstart = 298.15 "Start temperature of medium" annotation(Dialog(tab="General"));

      SI.Pressure pstart = 1e5 "Start pressure of medium" annotation(Dialog(tab="General"));

      Boolean useLiquidMassFraction = false "check if mass fraction for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMoleFraction and not useLiquidMolality),tab="General",group = "Reference composition"));

      Boolean useLiquidMoleFraction = false "check if molality for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMassFraction and not useLiquidMolality),tab="General",group = "Reference composition"));

      Boolean useLiquidMolality = false "check if molality for liquid species" annotation(choices(checkBox = true),Dialog(enable = (not useLiquidMoleFraction and not useLiquidMassFraction),tab="General",group = "Reference composition"));

      SI.MassFraction[nL] refXl = fill(1/nL,nL) "Reference mass fraction of liquid species (solutes+solvent) if useLiquidMassFraction" annotation(Dialog(enable = useLiquidMassFraction, tab="General",group = "Reference composition"));

      SI.MassFraction[nL] refYl = fill(1/nL,nL) "Reference mole fraction of liquid species (solutes+solvent) if useLiquidMoleFraction" annotation(Dialog(enable = useLiquidMoleFraction, tab="General",group = "Reference composition"));

      SI.MassFraction[nL-1] refml = fill(1e-5,nL-1) "Reference molality of liquid species (solutes) if useLiquidMolality" annotation(Dialog(enable = useLiquidMolality, tab="General",group = "Reference composition"));

      Boolean useSolidMassFraction = false "check if mass fraction for solid species" annotation(choices(checkBox = true),Dialog(enable = (not useSolidMoleFraction and not useSolidVolumeFraction),tab="General",group = "Reference composition"));

      Boolean useSolidMoleFraction = false "check if molality for solid species" annotation(choices(checkBox = true),Dialog(enable = (not useSolidMassFraction and not useSolidVolumeFraction),tab="General",group = "Reference composition"));

      Boolean useSolidVolumeFraction = false "check if mass fraction for solid species" annotation(choices(checkBox = true),Dialog(enable = (not useSolidMoleFraction and not useSolidMassFraction),tab="General",group = "Reference composition"));

      SI.MassFraction[ns] refXs = fill(1/ns,ns) "Reference mass fraction of solid species if useSolidMassFraction" annotation(Dialog(enable = useSolidMassFraction, tab="General",group = "Reference composition"));

      SI.MassFraction[ns] refYs = fill(1/ns,ns) "Reference mole fraction of solid species if useSolidMoleFraction" annotation(Dialog(enable = useSolidMoleFraction, tab="General",group = "Reference composition"));

      SI.MassFraction[ns] refPhis = fill(1/ns,ns) "Reference volume fraction of solid species if useSolidMoleFraction" annotation(Dialog(enable = useSolidVolumeFraction,  tab="General",group = "Reference composition"));

      Boolean usePhaseMassFraction = false "Check if liquid phase mass fraction" annotation(choices(checkBox = true),Dialog(enable = (not usePhaseMoleFraction and not usePhaseVolumeFraction),tab="General",group = "Reference composition"));

      Boolean usePhaseMoleFraction = false "Check if liquid phase mole fraction" annotation(choices(checkBox = true),Dialog(enable = (not usePhaseMassFraction and not usePhaseVolumeFraction),tab="General",group = "Reference composition"));

      Boolean usePhaseVolumeFraction = false "Check if liquid phase volume fraction" annotation(choices(checkBox = true),Dialog(enable = (not usePhaseMoleFraction and not usePhaseMassFraction),tab="General",group = "Reference composition"));

      Real refLiquidMassFraction = 0.5 "Reference liquid phase mass fraction" annotation(Dialog(enable = usePhaseMassFraction,tab="General",group = "Reference composition"));

      Real refLiquidMoleFraction = 0.5 "Reference liquid phase mole fraction" annotation(Dialog(enable = usePhaseMoleFraction,tab="General",group = "Reference composition"));

      Real refLiquidVolumeFraction = 0.5 "Reference liquid phase mass fraction" annotation(Dialog(enable = usePhaseVolumeFraction,tab="General",group = "Reference composition"));

      final SI.MassFraction[nL] Xl = if useLiquidMassFraction then refXl elseif useLiquidMoleFraction then moleToMassFractions(refYl,MMX[1+ns:ns+nL]) elseif useLiquidMolality then Functions.LiquidFunctions.calc_Xfromm(refml) else fill(1/nL,nL);

      final SI.MassFraction[ns] Xs = if useSolidMassFraction then refXs elseif useSolidMoleFraction then moleToMassFractions(refYs,MMX[1:ns]) elseif useSolidVolumeFraction then Functions.SolidFunctions.calc_XfromPhi(Tstart,pstart,refPhis) else fill(1/ns,ns);

      final SI.SpecificVolume vl = Functions.LiquidFunctions.calc_v(Tstart,pstart,Xl);

      final SI.SpecificVolume vs = Functions.SolidFunctions.calc_v(Tstart,pstart,Xs);

      final Real LiquidMassFractionFromVolume = refLiquidVolumeFraction*vs/((1-refLiquidVolumeFraction)*vl + refLiquidVolumeFraction*vs);

      final SI.MoleFraction[nL] Yl = massToMoleFractions(Xl,MMX[1+ns:ns+nL]);

      final SI.MoleFraction[ns] Ys = massToMoleFractions(Xs,MMX[1:ns]);

      final SI.MassFraction[ns+nL] refXfull = if usePhaseMassFraction then cat(1,(1-refLiquidMassFraction)*Xs,refLiquidMassFraction*Xl) elseif usePhaseMoleFraction then moleToMassFractions(cat(1,(1-refLiquidMoleFraction)*Ys,refLiquidMoleFraction*Yl),MMX)                                                                                                                                                                                                         elseif usePhaseVolumeFraction then cat(1,(1-LiquidMassFractionFromVolume)*Xs,LiquidMassFractionFromVolume*Xl) else fill(1/(ns+nL),ns+nL);

      final SI.MassFraction[nX] refX= calc_Xred(refXfull);

    end UserInterface;
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
    constant Real[nF,nX] lambda_id = Media.Common.Reaction.calc_lambda_id(nu_id) "lambda in identity transformation";// transpose(cat(2,-transpose(nu_id[:,nR+1:nF]),identity(nX)))
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

    function molarVolume "Calculates molar volume of mixture"
      input ThermodynamicState state;
      output MolarVolume v;
    protected
      Density d = density(state);
    algorithm
      v := molarMass(state)/d;
      annotation(Inline = true, smoothOrder = 3);
    end molarVolume;

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


    function T_phX "calculates temperature of gas"
      extends Modelica.Icons.Function;

      input SI.AbsolutePressure p "Pressure";
      input SI.SpecificEnthalpy h "Specific Enthalpy";
      input MassFraction[:] X;

      output SI.Temperature T "Temperature";

    protected
        package Internal
        "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
        extends Media.Common.OneNonLinearEquation;

         redeclare function extends f_nonlinear
         algorithm
         y :=Functions.calc_h(
                 x,
                 p_,
                 X_);
         end f_nonlinear;

        // Dummy definition has to be added for current Dymola
          redeclare function extends solve
          end solve;
        end Internal;

    protected
      SI.Temperature Tmax = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p);

    algorithm
      T := Internal.solve(
        h,
        273.15,
        Tmax,
        p,
        0,
        0,
        X);

      annotation(smoothOrder(normallyConstant = data)=20);
    end T_phX;

    function T_psX "calculates temperature of gas"
      extends Modelica.Icons.Function;

      input SI.AbsolutePressure p "Pressure";
      input SI.SpecificEntropy s "Specific Enthalpy";
      input MassFraction[:] X;

      output SI.Temperature T "Temperature";

    protected
        package Internal
        "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
        extends Media.Common.OneNonLinearEquation;

         redeclare function extends f_nonlinear
         algorithm
         y :=Functions.calc_s(
                 x,
                 p_,
                 X_);
         end f_nonlinear;

        // Dummy definition has to be added for current Dymola
          redeclare function extends solve
          end solve;
        end Internal;

    protected
      SI.Temperature Tmax = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p);

    algorithm
      T := Internal.solve(
        s,
        273.15,
        Tmax,
        p,
        0,
        0,
        X);

      annotation(smoothOrder(normallyConstant = data)=20);
    end T_psX;

    function p_TdX "calculates pressure of liquid mixture"
      extends Modelica.Icons.Function;

      input SI.Temperature T "Temperature";
      input SI.Density d "Pressure";
      input MassFraction[:] X;

      output SI.Pressure p;

    protected
        package Internal
        "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
        extends Media.Common.OneNonLinearEquation;

         redeclare function extends f_nonlinear
         algorithm
          y :=Functions.calc_d(
                 T_,
                 x,
                 X_);
         end f_nonlinear;

        // Dummy definition has to be added for current Dymola
          redeclare function extends solve
          end solve;
        end Internal;

    protected
        SI.Pressure pmin = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.psat(T);

         SI.SpecificVolume v = 1/d;

    algorithm
      p := Internal.solve(
        d,
        pmin,
        100e6,
        0,
        T,
        0,
        X);   //pmin,

      annotation(smoothOrder(normallyConstant = data)=20, experiment(Tolerance=1e-08));
    end p_TdX;

    function X_pTXred
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-20;
      parameter Real eps = 1e-8;

      MassFraction[nX] Xrednonzero;
      MassFraction[nX] Xred_;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      Integer index;

      Real[nX] Xs;
      Real[nF,nX] lambda_mass_id = P_to_id * lambda_mass;
    algorithm
      //ensure no species amount to be zero for numerical issues
      for i in 1:nX loop
        Xrednonzero[i] :=max(Modelica.Constants.eps, Xred[i]);
      end for;

      //convert reduced mass fraction basis
      Xs :=Modelica.Math.Matrices.solve(transpose(lambda_mass_id[nR+1:nF, :]), Xrednonzero);
      Xfull :=cat(1,zeros(nF - nX),Xs);
      Xfull :=P_to_orig*Xfull;

      Xred_ :=transpose(lambda_mass)*Xfull;

      // find initial guess x1 for Newton solver
      if sum(Xinit)> 0 then
      // ninit from initial guess provided as input (Xinit)
         scale :=sum(transpose(lambda_mass)*Xinit);
         minit :=Xinit*scale;
         ninit :=minit./MMX;
         for i in 1:nF loop
           x1[i] :=max(init, ninit[i]);
         end for;
      else
      // else generic initial guess
        for i in 1:nF loop
          x1[i] :=max(init, Xfull[i]/MMX[i]);
        end for;
      end if;

      // Newton algorithm
      solutionfound :=false;
      x1 :=P_to_id*x1;
      while not solutionfound loop
        f :=Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);
        x2 := Initialization.NewtonStep(x1,f,J,nF);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
          solutionfound:=true;
        end if;
        x1 :=x2;
      end while;
      x2 :=P_to_orig*x2;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      x_ :=Xfull;

    end X_pTXred;

    function XT_phXred
      "Liquid dissociation equilibrium calculated via Newton method with homotopy continuation method of step-by step adding solutes to initially pure aqueous solution"
      input Modelica.SIunits.Pressure p;
      input SI.SpecificEnthalpy h;
      input MassFraction[nX] Xred;
      output Real[nF] X;
      output Real T;

    protected
      parameter Real eps = 1e-5;

      SI.Temperature Tmax = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p);
      SI.Temperature Tmin = 273.15;

      Real[1] f;
      Real[1,1] J;
      SI.Temperature T_;
      SI.Temperature[1] Ttemp;

    algorithm
      //calculate start values for inner-outer algorithm
      T :=exp((log(Tmin) + log(Tmax))/2);
      X :=X_pTXred(p,T,Xred);

      //inner (equilibrium) - outer (enthalpy) algorithm
      while Modelica.Math.Vectors.norm({T_/T-1}) > eps loop
        f := Initialization.calc_f.Enthalpy(T,p,h,X);
        J := Initialization.calc_J.Enthalpy(T,p,X);
        Ttemp:=Initialization.NewtonStep({T},f,J,1);
        assert(sum(Ttemp)>0,"Error in Newton step for enthalpy calculation",AssertionLevel.error);
        T_ :=T;
        T :=Ttemp[1];

        if T <Tmin then
          T :=Tmin;
          T_ :=exp((log(Tmin) + log(Tmax))/2);
        end if;

        //inner algorithm
        X :=X_pTXred(p,T,Xred,X);
    end while;

    end XT_phXred;

    function XT_psXred
      "Liquid dissociation equilibrium calculated via Newton method with homotopy continuation method of step-by step adding solutes to initially pure aqueous solution"
      input Modelica.SIunits.Pressure p;
      input SI.SpecificEntropy s;
      input MassFraction[nX] Xred;
      output Real[nF] X;
      output Real T;

    protected
      parameter Real eps = 1e-5;

      SI.Temperature Tmax = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p);
      SI.Temperature Tmin = 273.15;

      Real[1] f;
      Real[1,1] J;
      SI.Temperature T_;
      SI.Temperature[1] Ttemp;

    algorithm
      //calculate start values for inner-outer algorithm
      T :=exp((log(Tmin) + log(Tmax))/2);
      X :=X_pTXred(p,T,Xred);

      //inner (equilibrium) - outer (enthalpy) algorithm
      while Modelica.Math.Vectors.norm({T_/T-1}) > eps loop
        f := Initialization.calc_f.Entropy(T,p,s,X);
        J := Initialization.calc_J.Entropy(T,p,X);
        Ttemp:=Initialization.NewtonStep({T},f,J,1);
        assert(sum(Ttemp)>0,"Error in Newton step for enthalpy calculation",AssertionLevel.error);
        T_ :=T;
        T :=Ttemp[1];

        if T <Tmin then
          T :=Tmin;
          T_ :=exp((log(Tmin) + log(Tmax))/2);
        end if;

        //inner algorithm
        X :=X_pTXred(p,T,Xred,X);
      end while;

    end XT_psXred;

    function Xp_dTXred
      "Gas-liquid and liquid dissociation equilibrium calculated via Newton method with homotopy continuation method of step-by step adding solutes to initially pure aqueous solution"
      input Modelica.SIunits.Density d;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      output MassFraction[ns+nL] X;
      output SI.Pressure p;

    protected
      parameter Real eps = 1e-5;

      SI.Pressure pmin = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.psat(T);
      SI.Pressure pmax = 1000e5;

      Real[1] f={1};
      Real[1,1] J;
      SI.Pressure p_=pmin;
      SI.Pressure[1] ptemp;

      SI.SpecificVolume v = 1/d;

    algorithm
      //calculate start values for inner-outer algorithm
      p :=exp((log(pmin) + log(pmax))/2);//1.1*pmin;//
      X :=X_pTXred(p,T,Xred);

      //inner (equilibrium) - outer (enthalpy) algorithm
      while Modelica.Math.Vectors.norm({p_/p-1},Modelica.Constants.inf) > eps loop

        f := Initialization.calc_f.Density(T,p,X,d);
        J := Initialization.calc_J.Density(T,p,X);
        ptemp:=Initialization.NewtonStep({p},f,J,1);
        assert(sum(ptemp)>0,"Error in Newton step for enthalpy calculation",AssertionLevel.error);
        p_ :=p;
        p :=ptemp[1];

       if p <pmin then
         p :=pmin;
         p_ :=exp((log(pmin) + log(pmax))/2);
       end if;

        //inner algorithm
        X :=X_pTXred(p,T,Xred,X);
      end while;
    end Xp_dTXred;

    function calc_Xred
      input SI.MassFraction[ns+nL] Xfull;
      output SI.MassFraction[nX] Xred;

    algorithm
      Xred :=transpose(lambda_mass)*Xfull/(sum(transpose(lambda_mass)*Xfull));

    end calc_Xred;

    package Initialization "Package containing all functions regarding initialization"
    //   constant Integer ns = 3;
    //   constant DataRecordIG[ns] dataIG;
    //   constant DataRecordIG[ns] dataPR;
    //   constant Media.Common.Types.GasModel GasModel = MixtureLiquid.GasModel;
    //
    //   constant Integer nL=2;
    //   final constant Integer nLi = nL-1;
    //   constant DataRecordL[nL-1] datal;
    //   constant Media.Common.Types.LiquidModel LiquidModel;

      package calc_f

        function Enthalpy
          "Liquid dissociation equilibrium with constraints on solute molalities"
          input SI.Temperature T;
          input SI.Pressure p;
          input SI.SpecificEnthalpy h;
          input MassFraction[ns+nL] X;
          output Real[1] f;

        algorithm
          //Enthalpy balance
          f[1] :=(Functions.calc_h(T,p,X) - h)/Common.StdRefH2O.h_tr;
        end Enthalpy;

        function Entropy
          "Liquid dissociation equilibrium with constraints on solute molalities"
          input SI.Temperature T;
          input SI.Pressure p;
          input SI.SpecificEntropy s;
          input MassFraction[ns+nL] X;
          output Real[1] f;

        algorithm
          f[1] :=(Functions.calc_s(T,p,X) - s)/Common.StdRefH2O.s_tr;
        end Entropy;

        function Density
          "Liquid dissociation equilibrium with constraints on solute molalities"
          input SI.Temperature T;
          input SI.Pressure p;
          input MassFraction[ns+nL] X;
          input SI.Density d;
          output Real[1] f;

        algorithm
           f[1] :=(Functions.calc_d(T,p,X) - d);
        end Density;

        function SLE_Tp
          "SLE initialization according to Leal et al. (2016): calculation of residuals"
          input Real[nF] x;
          input Real[nX] Xred;
          input SI.Temperature T;
          input SI.Pressure p;
          output Real[nF] f;

        protected
          SI.MolarEnergy gi[ns+nL];
          SI.MolarEnergy gi_id[ns+nL];
          SI.MolarEnergy gr[nR];
          Real[nR] logK;
          Real[nL] gamma_i;
          Real[nL] m;
          Real [nL] a;
          SI.MoleFraction[ns+nL] Y;
          SI.MassFraction[ns+nL] X;
          SI.MoleFraction[nL] Yl;
          SI.MassFraction[nL] Xl;
          SI.MoleFraction[ns] Ys;
          SI.MassFraction[ns] Xs;
          Real tau = 1e-37;
          SI.MassFraction[ns+nL] z;// = {if x[i] > 0 then tau/max(tau,x[i]) else 0 for i in 1:ns+nL};//tau./x;

          Real[ns+nL] logreacBase;
          Real[ns+nL] logreacBase_id;
          Real[ns+nL] x_orig;

          SI.Mass[ns+nL] mass;

        algorithm

          //assert(sum(x[1:ns+nL]) > 0, "x is zero in SLE_Tp");

          x_orig :=P_to_orig*x;
          for i in 1:ns+nL loop
            x_orig[i] :=x_orig[i];//max(tau, x_orig[i]);//x_orig without pivoting, max formulation with pivoting
          end for;
        //     z :={if x_orig[i] > 0 then tau/max(1e-20*tau, x_orig[i]) else 0 for i in 1:ns + nL};
               z :=tau ./ x_orig;

          Y :=x_orig[1:ns+nL]/sum(x_orig[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          Xs :=Functions.calc_X_M(Ys, MMX[1:ns]);
          Xl :=Functions.calc_X_M(Yl, MMX[1 + ns:ns + nL]);

          m[1:nL-1] := Functions.LiquidFunctions.calc_mfromY(Yl);
          m[nL] := Yl[nL]/MH2O;

          gamma_i := Functions.LiquidFunctions.calc_gamma(T,p,Xl);
          gi := Functions.calc_gi(T, p) .* MMX;
          gi_id :=P_to_id*gi;

        //    gr :=nu_id*gi_id;
           gr :=nu*gi;
          logK :=-gr/Modelica.Constants.R/T;

          a :=gamma_i .* m;
          logreacBase :=cat(1,-z[1:ns],log(a)-z[ns+1:ns+nL]);
           logreacBase_id :=P_to_id*logreacBase;

        //   x_id :=P_to_id*x;

          //isopotential
           f[1:nR] :=nu*logreacBase - logK;
        //   f[1:nR] :=nu_id*logreacBase_id - logK;

          //reduced mass balance
          // f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
        //   f[nR+1:nF] :=transpose(lambda_nu)*x_id - Xred;//transpose(lambda_mass)*mass - Xred;
          f[nR+1:nF] :=transpose(lambda_id)*x - Xred;//transpose(lambda_mass)*mass - Xred;
        end SLE_Tp;

      end calc_f;

      package calc_J

        function Enthalpy
          "Liquid dissociation equilibrium with constraints on solute molalities"
          input SI.Temperature T;
          input SI.Pressure p;
          input MassFraction[ns+nL] X;
          output Real[1,1] J;

        protected
          SI.SpecificHeatCapacity cp= Functions.calc_cp(T,p, X);
        algorithm
          J[1,1] :=cp/Common.StdRefH2O.h_tr;///h

        end Enthalpy;

        function Entropy
          "Liquid dissociation equilibrium with constraints on solute molalities"
          input Temperature T;
          input SI.Pressure p;
          input MassFraction[ns+nL] X;
          output Real[1,1] J;

        protected
          SI.SpecificHeatCapacity cp=Functions.calc_cp(T,p,X);
        algorithm
          J[1,1] :=cp/T/Common.StdRefH2O.s_tr;

        end Entropy;

        function Density
          "Liquid dissociation equilibrium with constraints on solute molalities"
          input SI.Temperature T;
          input SI.Pressure p;
          input MassFraction[ns+nL] X;
          output Real[1,1] J;
        protected
          SI.DerDensityByPressure ddpT=Functions.calc_ddpT(T,p,X);
          SI.Density d=Functions.calc_d(T,p,X);
        algorithm
          J[1,1] :=ddpT;

        end Density;

        function SLE_Tp
          "SLE initialization according to Leal et al. (2016): calculation of Jacobian"
          input Real[nF] x;
          input SI.Temperature T;
          input SI.Pressure p;
          output Real[nF,nF] J;

        protected
          SI.AmountOfSubstance[ns+nL] n;
          SI.MoleFraction[ns+nL] Y;
          SI.MassFraction[ns+nL] X;
          SI.MoleFraction[nL] Yl;
          SI.MassFraction[nL] Xl;
          SI.MoleFraction[ns] Ys;
          SI.MassFraction[ns] Xs;
          Real tau = 1e-37;
          SI.MassFraction[ns+nL] z;
          Real[ns+nL,ns+nL] H;
          Real[ns+nL,ns+nL] H_id;
          Real[ns+nL] diagH;
          Real[ns+nL] x_orig;
        algorithm

          x_orig :=P_to_orig*x;

          for i in 1:ns+nL loop
            n[i] :=x_orig[i];//max(tau, x_orig[i]);//x_orig without pivoting, max formulation with pivoting
          end for;

          Y :=n[1:ns+nL]/sum(n[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          z :=tau./n;

          diagH[1:ns] :=tau ./ n[1:ns] .^ 2;//z[1:ns]./n[1:ns];//
          diagH[ns+1:ns+nL] :=(1 .- Yl .+ z[ns+1:ns+nL]) ./ n[ns + 1:ns + nL];// + tau ./ n[ns+1:ns+nL] .^ 2;//(1 .- Yl) ./ n[ns + 1:ns + nL] + tau ./ n[ns+1:ns+nL] .^ 2;

          H :=diagonal(diagH);//fill(diagH,ns+nL);//
          H_id :=transpose(P_to_id)*H;

        //    for i in 1:ns+nL loop
        //      //isopotential
        //      for r in 1:nR loop
        //        if i < ns+1 then
        //           J[r,i] := nu[r,i]*tau/x[i]^2 else 0;//nu[r,i]*z[i]/n[i] else 0;//if abs(nu[r,i]) > 0 then nu[r,i]*z[i]/n[i] else 0;
        //         else
        //           J[r,i] := if abs(nu[r,i])>0 then nu[r, i]*(1 - Yl[i - ns])/x[i] + nu[r,i]*tau/x[i]^2 else 0;//nu[r, i]*(1 - Yl[i - ns])/n[i] + nu[r,i]*z[i]/n[i] else 0;//if abs(nu[r,i]) > 0 then nu[r, i]*(1 - Yl[i - ns])/n[i] + nu[r,i]*z[i]/n[i] else 0;
        //         end if;
        //      end for;
        //    end for;

          //isopotential
        //   for r in 1:nR loop
        //     J[r,:] :=nu[r,:].*diagH;
        //   end for;
             J[1:nR,:] :=nu*H;//nu_id__ .* H;//nu_id__*H;//
             J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id);//transpose(P_nu_id);

        //      J[1:nR,:] :=nu_id_ .* H_id;
        //     J[1:nR,:] :=nu_id .* H_id;

          //reduced mass balance
           J[nR+1:nF,:] :=transpose(lambda_id);
        end SLE_Tp;

      end calc_J;

      function NewtonStep
        "Calculates Newton step with step size such that x[t,i] > 0"
        input Real[nNS] x;
        input Real[nNS] f;
        input Real[nNS,nNS] J;
        input Integer nNS;
        output Real[nNS] x_;

      protected
        Real beta=1;
        Real[nNS] alpha;
        Real delta = 0.9999;
         Real[nNS,nNS] invJ= Modelica.Math.Matrices.inv(J);
        Real[nNS] delta_x =  invJ*f;//Modelica.Math.Matrices.solve(J,f);//
      algorithm

         for i in 1:nNS loop
           if x[i] - delta_x[i] < 0 then
             alpha[i] :=delta*x[i]/delta_x[i];
           else
             alpha[i] :=1;
           end if;
         end for;
         beta :=min(alpha);
         x_ :=x - beta*delta_x;

      //    while min(x_temp) < 0 loop
      //      x_temp :=x - 1/beta*invJ*f;
      //      beta :=beta *10;
      //    end while;
      //    assert(min(x_temp)>0,"Newton step too small");
      //    x_ :=x_temp;

      end NewtonStep;
    end Initialization;

    package Functions
      extends Common.Functions(nSfunfun= ns, dataSfunfun = datas, nLfunfun = nL, dataLfunfun = datal, LiquidModelfunfun = LiquidModel, interactionLfunfun = interactionL);
    end Functions;

  end MixtureLiquid;

end Common;
