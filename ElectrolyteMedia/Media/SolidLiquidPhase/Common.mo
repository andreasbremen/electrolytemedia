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
  //      nXi = nX-1,
  //      fixedX=false,
  //      reducedX = true,
  //      singleState=false,

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

      Real[nR,ns+nL] nu = zeros(nR,ns+nL) "Stoichiometry matrix of solid-liquid and dissociation equilibria" annotation(Dialog(tab="Reaction"));

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
    constant Real[nR,nF] nu_id_= Media.Common.Reaction.calc_nu_id(nu);
                                                                      //nu
    constant Integer[nF,nF] P_to_orig_=Media.Common.Reaction.calc_P_nu_id(nu);
                                                                              //identity(nF);
    constant Integer[nF,nF] P_to_id_ = transpose(P_to_orig_);
    constant Real[nR,nF] nu_id__ = nu_id_*transpose(P_to_orig_) "nu in original order";
    constant Real[nF,nX] lambda_nu = transpose(cat(2,-transpose(nu_id_[:,nR+1:nF]),identity(nX)));
  //   constant Real[nF,nX] lambda_id = cat(2,identity(nX),transpose(nu_id[1:nX,nX+1:nF]));
    constant Real[nR,nF] nu_mass = Functions.Reaction.calc_nu_mass(nu);
    constant Real[nF] MMX_ = P_to_id_*MMX;

    constant UserInterface
      userInterface "User interface"
      annotation (Placement(transformation(extent={{-10,-8},{10,12}})));

    constant Real[nF,nX] lambda_=
      Media.Common.Reaction.calc_lambda(
      nu,
      nR,
      nF);
    constant Real[nF,nX] lambda_id = Media.Common.Reaction.calc_lambda_id(lambda_);
    constant Real[nF,nX] lambda = P_to_orig*lambda_id;
    constant Integer[nF,nF] P_to_orig= Media.Common.Reaction.calc_P_lambda_id(lambda_);
    constant Integer[nF,nF] P_to_id = transpose(P_to_orig);
    constant Real[nR,nF] nu_id = nu*transpose(P_to_id);
    constant Real[nF,nX] lambda_mass_=P_to_orig_*{{lambda_nu[i,j]/MMX_[i] for j in 1:nX} for i in 1:nF} "Null space of mass based stoichiometry matrix";
    constant Real[nF,nX] lambda_mass={{lambda_[i,j]/MMX[i] for j in 1:nX} for i in 1:nF} "Null space of mass based stoichiometry matrix";
  //  constant Real[nX] scale = {sum(lambda_mass_[i,j] for i in 1:nF) for j in 1:nX};
  //   constant Real[:,:] lambda_mass = {{lambda_mass_[i,j]/scale[j] for j in 1:nX} for i in 1:nF};

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
      parameter Real[nX] Xredstart_ = transpose(lambda_mass_)*Xfullstart "for scaling of species moles consistent with maxx invariant vector X";
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
      Real[nF,nX] lambda_mass_id = P_to_id_ * lambda_mass;
    algorithm
      //ensure no species amount to be zero for numerical issues
      for i in 1:nX loop
        Xrednonzero[i] :=max(Modelica.Constants.eps, Xred[i]);
      end for;

      //convert reduced mass fraction basis
      Xs :=Modelica.Math.Matrices.solve(transpose(lambda_mass_id[nR+1:nF, :]), Xrednonzero);
      Xfull :=cat(1,zeros(nF - nX),Xs);
      Xfull :=P_to_orig_*Xfull;

      Xred_ :=transpose(lambda_mass_)*Xfull;

      // find initial guess x1 for Newton solver
      if sum(Xinit)> 0 then
      // ninit from initial guess provided as input (Xinit)
         scale :=sum(transpose(lambda_mass_)*Xinit);
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
      x1 :=P_to_id_*x1;
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
      x2 :=P_to_orig_*x2;
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

          x_orig :=P_to_orig_*x;
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
          gi_id :=P_to_id_*gi;

        //    gr :=nu_id*gi_id;
           gr :=nu_id__*gi;
          logK :=-gr/Modelica.Constants.R/T;

          a :=gamma_i .* m;
          logreacBase :=cat(1,-z[1:ns],log(a)-z[ns+1:ns+nL]);
           logreacBase_id :=P_to_id_*logreacBase;

        //   x_id :=P_to_id_*x;

          //isopotential
           f[1:nR] :=nu_id__*logreacBase - logK;
        //   f[1:nR] :=nu_id*logreacBase_id - logK;

          //reduced mass balance
          // f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
        //   f[nR+1:nF] :=transpose(lambda_nu)*x_id - Xred;//transpose(lambda_mass)*mass - Xred;
          f[nR+1:nF] :=transpose(lambda_nu)*x - Xred;//transpose(lambda_mass)*mass - Xred;
        end SLE_Tp;

        function SLE_Tp_
          "SLE initialization according to Leal et al. (2016): calculation of residuals"
          input Real[nF] x;
          input Real[nX] Xred;
          input SI.Temperature T;
          input SI.Pressure p;
          output Real[nF] f;

        protected
          SI.MolarEnergy gi[ns+nL];
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
          Real[ns] slack = 1e-20./x[1:ns];//x[nF+1:nF+ns];
          parameter Real eps = 1e-15;

          Real[ns+nL] logreacBase;

          SI.Mass[ns+nL] mass;

        algorithm
          //assert(sum(x[1:ns+nL]) > 0, "x is zero in SLE_Tp");

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          Xs :=Functions.calc_X_M(Ys, MMX[1:ns]);
          Xl :=Functions.calc_X_M(Yl, MMX[1 + ns:ns + nL]);

          m[1:nL-1] := Functions.LiquidFunctions.calc_mfromY(Yl);
          m[nL] := Yl[nL]/MH2O;

          gamma_i := Functions.LiquidFunctions.calc_gamma(T,p,Xl);
          gi := Functions.calc_gi(T, p) .* MMX;

          gr :=nu*gi;
          logK :=-gr/Modelica.Constants.R/T;

          a :=gamma_i .* m;
          logreacBase :=cat(1,-slack,log(a));

          mass :=x[1:nF] .* MMX;

          //isopotential
          f[1:nR] :=nu*logreacBase - logK;

          //reduced mass balance
          f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
        end SLE_Tp_;

        function SLE_Tp_v1
          "SLE initialization according to Leal et al. (2016): calculation of residuals"
          input Real[nF] x;
          input Real[nX] Xred;
          input SI.Temperature T;
          input SI.Pressure p;
          output Real[nF] f;

        protected
          SI.MolarEnergy gi[ns+nL];
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
          Real tau = 1e-30;
          SI.MassFraction[ns+nL] z = {if x[i] > 0 then tau/max(tau,x[i]) else 0 for i in 1:ns+nL};//tau./x;

          Real[ns+nL] logreacBase;

          SI.Mass[ns+nL] mass;

        algorithm
          assert(sum(x[1:ns+nL]) > 0, "x is zero in SLE_Tp");

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          Xs :=Functions.calc_X_M(Ys, MMX[1:ns]);
          Xl :=Functions.calc_X_M(Yl, MMX[1 + ns:ns + nL]);

          m[1:nL-1] := Functions.LiquidFunctions.calc_mfromY(Yl);
          m[nL] := Yl[nL]/MH2O;

          gamma_i := Functions.LiquidFunctions.calc_gamma(T,p,Xl);
          gi := Functions.calc_gi(T, p) .* MMX;

          gr :=nu*gi;
          logK :=-gr/Modelica.Constants.R/T;

          a :=gamma_i .* m;
          logreacBase :=cat(1,-z[1:ns],log(a)-z[ns+1:ns+nL]);

          mass :=x[1:nF] .* MMX;

          //isopotential
          f[1:nR] :=nu*logreacBase - logK;

          //reduced mass balance
          // f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
          f[nR+1:nF] :=transpose(lambda_)*x - Xred;//transpose(lambda_mass)*mass - Xred;
        end SLE_Tp_v1;

        function SLE_Tp_v2
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
          Real tau = 1e-30;
          SI.MassFraction[ns+nL] z = {if x[i] > 0 then tau/max(tau,x[i]) else 0 for i in 1:ns+nL};//tau./x;

          Real[ns+nL] logreacBase;
          Real[ns+nL] logreacBase_id;
          Real[ns+nL] x_id;

          SI.Mass[ns+nL] mass;

        algorithm
          assert(sum(x[1:ns+nL]) > 0, "x is zero in SLE_Tp");

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
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

          mass :=x[1:nF] .* MMX;
          x_id :=P_to_id*x;

          //isopotential
           f[1:nR] :=nu*logreacBase - logK;
        //   f[1:nR] :=nu_id*logreacBase_id - logK;

          //reduced mass balance
          // f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
          f[nR+1:nF] :=transpose(lambda_id)*x_id - Xred;//transpose(lambda_mass)*mass - Xred;
        end SLE_Tp_v2;

        function SLE_Tp_v3
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
          Real tau = 1e-30;
          SI.MassFraction[ns+nL] z;// = {if x[i] > 0 then tau/max(tau,x[i]) else 0 for i in 1:ns+nL};//tau./x;

          Real[ns+nL] logreacBase;
          Real[ns+nL] logreacBase_id;
          Real[ns+nL] x_orig;

          SI.Mass[ns+nL] mass;

        algorithm

          assert(sum(x[1:ns+nL]) > 0, "x is zero in SLE_Tp");

          x_orig :=P_to_orig_*x;

          z :={if x_orig[i] > 0 then tau/max(tau, x_orig[i]) else 0 for i in 1:ns + nL};

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
          gi_id :=P_to_id_*gi;

        //    gr :=nu_id*gi_id;
           gr :=nu_id__*gi;
          logK :=-gr/Modelica.Constants.R/T;

          a :=gamma_i .* m;
          logreacBase :=cat(1,-z[1:ns],log(a)-z[ns+1:ns+nL]);
           logreacBase_id :=P_to_id_*logreacBase;

        //   x_id :=P_to_id_*x;

          //isopotential
           f[1:nR] :=nu_id__*logreacBase - logK;
        //   f[1:nR] :=nu_id*logreacBase_id - logK;

          //reduced mass balance
          // f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
        //   f[nR+1:nF] :=transpose(lambda_nu)*x_id - Xred;//transpose(lambda_mass)*mass - Xred;
          f[nR+1:nF] :=transpose(lambda_nu)*x - Xred;//transpose(lambda_mass)*mass - Xred;
        end SLE_Tp_v3;
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

          x_orig :=P_to_orig_*x;

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
          H_id :=transpose(P_to_id_)*H;

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
             J[1:nR,:] :=nu_id__*H;//nu_id__ .* H;//nu_id__*H;//
             J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id_);//transpose(P_nu_id);

        //      J[1:nR,:] :=nu_id_ .* H_id;
        //     J[1:nR,:] :=nu_id .* H_id;

          //reduced mass balance
           J[nR+1:nF,:] :=transpose(lambda_nu);
        end SLE_Tp;

        function SLE_Tp_
          "SLE initialization according to Leal et al. (2016): calculation of Jacobian"
          input Real[nF] x;
          input SI.Temperature T;
          input SI.Pressure p;
          output Real[nF,nF] J;

        protected
          SI.MolarEnergy gi[ns+nL];
          SI.MoleFraction[ns+nL] Y;
          SI.MassFraction[ns+nL] X;
          SI.MoleFraction[nL] Yl;
          SI.MassFraction[nL] Xl;
          SI.MoleFraction[ns] Ys;
          SI.MassFraction[ns] Xs;
          parameter Real eps = 1e-15;
          Real[ns] slack = 1e-20./x[1:ns];//x[nF+1:nF+ns];
          Real[ns+nL] logreacBase;
        algorithm

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          for i in 1:ns+nL loop
            //isopotential
            for r in 1:nR loop
              if i < ns+1 then
                J[r,i] := if x[i] > 0 then nu[r,i]*slack[i]/x[i] else 0;
               else
                 J[r,i] :=if nu[r, i] > 0 then nu[r, i]*(1 - Yl[i - ns])/x[i] else 0;
               end if;
            end for;

            //reduced mass balance
            for j in 1:ns+nL - nR loop
              J[j + nR,i] := lambda_mass[i, j]*MMX[i];
            end for;
          end for;
        end SLE_Tp_;

        function SLE_Tp_v1
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
          Real tau = 1e-20;
          SI.MassFraction[ns+nL] z;
          Real[ns+nL,ns+nL] H;
          Real[ns+nL] diagH;
        algorithm

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          for i in 1:ns+nL loop
            n[i] :=max(tau, x[i]);
          end for;
          z :=tau./n;

          diagH[1:ns] :=tau ./ n[1:ns] .^ 2;
          diagH[ns+1:ns+nL] :=(1 .- Yl) ./ n[ns + 1:ns + nL] + tau ./ n[ns+1:ns+nL] .^ 2;

          H :=fill(diagH,ns+nL);

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
          J[1:nR,:] :=nu .* H;

          //reduced mass balance
           J[nR+1:nF,:] :=transpose(lambda_);
        end SLE_Tp_v1;

        function SLE_Tp_v2
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
          Real tau = 1e-20;
          SI.MassFraction[ns+nL] z;
          Real[ns+nL,ns+nL] H;
          Real[ns+nL,ns+nL] H_id;
          Real[ns+nL] diagH;
        algorithm

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          for i in 1:ns+nL loop
            n[i] :=max(tau, x[i]);
          end for;
          z :=tau./n;

          diagH[1:ns] :=tau ./ n[1:ns] .^ 2;
          diagH[ns+1:ns+nL] :=(1 .- Yl) ./ n[ns + 1:ns + nL] + tau ./ n[ns+1:ns+nL] .^ 2;

          H :=fill(diagH,ns+nL);
          H_id :=P_to_id*H;

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
             J[1:nR,:] :=nu .* H;
             J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id);
        //     J[1:nR,:] :=nu_id .* H_id;

          //reduced mass balance
           J[nR+1:nF,:] :=transpose(lambda_id);
        end SLE_Tp_v2;

        function SLE_Tp_v3
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
          Real tau = 1e-20;
          SI.MassFraction[ns+nL] z;
          Real[ns+nL,ns+nL] H;
          Real[ns+nL,ns+nL] H_id;
          Real[ns+nL] diagH;
          Real[ns+nL] x_orig;
        algorithm

          x_orig :=P_to_orig_*x;

          Y :=x_orig[1:ns+nL]/sum(x_orig[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          for i in 1:ns+nL loop
            n[i] :=max(tau, x_orig[i]);
          end for;
          z :=tau./n;

          diagH[1:ns] :=tau ./ n[1:ns] .^ 2;
          diagH[ns+1:ns+nL] :=(1 .- Yl) ./ n[ns + 1:ns + nL] + tau ./ n[ns+1:ns+nL] .^ 2;

          H :=fill(diagH,ns+nL);
          H_id :=P_to_id_*H;

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
             J[1:nR,:] :=nu_id__ .* H;
             J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id_);//transpose(P_nu_id);
        //     J[1:nR,:] :=nu_id .* H_id;

          //reduced mass balance
           J[nR+1:nF,:] :=transpose(lambda_nu);
        end SLE_Tp_v3;
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

    function X_pTXred_
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-8;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B;
      Real[nR] b;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      Integer index;
    algorithm
        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B[r,firstnonzero[r]] :=1;
           b[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xred,B,b);
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;

      // Newton algorithm
      solutionfound :=false;
      while not solutionfound loop
        f :=Initialization.calc_f.SLE_Tp(x1,Xred,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);
        x2 := Initialization.NewtonStep(x1,f,J,nF);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
          solutionfound:=true;
        end if;
        x1 :=x2;
      end while;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      for i in 1:nF loop
      x_[i] :=max(1e-20,Xfull[i]);
      end for;

    end X_pTXred_;

    function X_pTXred_v2
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-8;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B;
      Real[nR] b;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      MassFraction[nX] Xred_;

      Integer index;
    algorithm
        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B[r,firstnonzero[r]] :=1;
           b[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xred,B,b);
         assert(sum(Xfull_in)> 0,"Test");
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;

         Xred_ :=transpose(lambda_mass_)*Xfull_in;
      // Newton algorithm
      solutionfound :=false;
      while not solutionfound loop
    //     f :=Initialization.calc_f.SLE_Tp_v1(x1,Xred,T,p);
    //     J :=Initialization.calc_J.SLE_Tp_v1(x1,T,p);
    //     x2 := Initialization.NewtonStep(x1,f,J,nF);

        f :=Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);
        x2 := P_to_orig_*Initialization.NewtonStep(P_to_id_*x1,f,J,nF);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
          solutionfound:=true;
        end if;
        x1 :=x2;
      end while;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      for i in 1:nF loop
      x_[i] :=max(1e-20,Xfull[i]);
      end for;

    end X_pTXred_v2;

    function X_pTXred_v3
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-8;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B_;
      Real[nR] b_;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      MassFraction[nX] Xred_;

      Integer index;
      Boolean success;
      Real[nF,nX] lambda_mass1;

      Real[nR,nR] D;
      Real[nR,nX] B;
      Real[nX,nR] C;
      Real[nX,nX] I;
      Real[nR] a;
      Real[nX] b;
      Real[nR,nR] invD;
      Integer nNS;
    algorithm
        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B_[r,firstnonzero[r]] :=1;
           b_[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xred,B_,b_);
         assert(sum(Xfull_in)> 0,"Test");
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;

         lambda_mass1 :=P_to_orig_*lambda_mass;
         // Xred_ :=Modelica.Math.Matrices.solve(transpose(lambda_mass),Xred);//
         Xred_:=transpose(lambda_mass_)*Xfull_in;
         success :=Modelica.Utilities.Streams.writeRealMatrix(
        "U:/OutputX_pTXred.mat",
        "J",
        lambda_mass1);
      // Newton algorithm
      solutionfound :=false;
      x1 :=P_to_id_*x1;
      while not solutionfound loop
    //     f :=Initialization.calc_f.SLE_Tp_v1(x1,Xred,T,p);
    //     J :=Initialization.calc_J.SLE_Tp_v1(x1,T,p);
    //     x2 := Initialization.NewtonStep(x1,f,J,nF);

        f :=Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);
        a :=f[1:nR];
        b :=f[nR + 1:nF];
        D :=J[1:nR, 1:nR];
        B :=J[1:nR, nR + 1:nF];
        C :=J[nR + 1:nF, 1:nR];
        I :=identity(nX);
        invD :={{if i == j then 1/D[i, j] else 0 for i in 1:nR} for j in 1:nR};
        x2[nR+1:nF] :=Initialization.NewtonStep(
          x1[nR + 1:nF],
          b - C*invD*a,
          I - C*invD*B,
          nX);
        x2[1:nR] :=invD*(a - B*x2[nR + 1:nF]);
    //     x2 := P_to_orig_*Initialization.NewtonStep(P_to_id_*x1,f,J,nF);
         x2 := Initialization.NewtonStep(x1,f,J,nF);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
          solutionfound:=true;
        end if;
        x1 :=x2;
      end while;
      x2 :=P_to_orig_*x2;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      for i in 1:nF loop
      x_[i] :=max(1e-20,Xfull[i]);
      end for;

    end X_pTXred_v3;

    function X_pTXred_v4
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-3;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;
      Boolean solutionfound1 = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B_;
      Real[nR] b_;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      MassFraction[nX] Xred_;

      Integer index;
      Boolean success;
      Real[nF,nX] lambda_mass1;

      Real[nR,nR] D;
      Real[nR,nR] Dp;
      Real[nR,nR] Dn;
      Real[nR,nR] invDp;
      Real[nR,nR] invDn;
      Real[nR] Dvector;
      Real[nR] Dpvector;
      Real[nR] Dnvector;
      Real[nR] invDvector;
      Real[nR,nX] B;
      Real[nR,nX] Bp;
      Real[nR,nX] Bn;
      Real[nX,nR] C;
      Real[nX,nX] I;
      Real[nR] a;
      Real[nX] b;
      Real[nR,nR] invD;
      Boolean[nR] I_p;
      Integer n_p;
      Integer n_n;
      Integer[nR] pivots;
      Integer[nR,nR] P;
      Real[nF,nF] J_;
      Real[nF] f_;

      Real[nR] ap;
      Real[nR] an;
      Real[nX,nR] Cp;
      Real[nX,nR] Cn;
    algorithm
        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B_[r,firstnonzero[r]] :=1;
           b_[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xred,B_,b_);
         assert(sum(Xfull_in)> 0,"Test");
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;
    //        x1[1:nF] :=cat(1,{0.001*i for i in 1:nF-1},{1-0.001*sum(i for i in 1:nF-1)});//fill(1/nF, nF);
        x1[1:nF] :=fill(1/nF, nF);

         lambda_mass1 :=P_to_orig_*lambda_mass;
         // Xred_ :=Modelica.Math.Matrices.solve(transpose(lambda_mass),Xred);//
         Xred_:=transpose(lambda_mass_)*Xfull_in;

      // Newton algorithm
      solutionfound :=false;
      x1 :=P_to_id_*x1;
       index :=1;
        f :=Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);
      while not solutionfound loop
    //     f :=Initialization.calc_f.SLE_Tp_v1(x1,Xred,T,p);
    //     J :=Initialization.calc_J.SLE_Tp_v1(x1,T,p);
    //     x2 := Initialization.NewtonStep(x1,f,J,nF);
    //     assert(sum(x1)>0, "x1 is zero at iteration " + String(index));

         a :=f[1:nR];
         b :=f[nR + 1:nF];
         D :=J[1:nR, 1:nR];
         B :=J[1:nR, nR + 1:nF];
         C :=J[nR + 1:nF, 1:nR];
         I :=identity(nX);
         I_p :={false for i in 1:nR};//{abs(D[i, i]) >= Modelica.Math.Vectors.norm(C[:, i],Modelica.Constants.inf) for i in 1:nR};

         n_p :=Modelica.Math.BooleanVectors.countTrue(I_p);
         pivots[1:n_p] :=Modelica.Math.BooleanVectors.index(I_p);
         pivots[n_p+1:nR] :=Modelica.Math.BooleanVectors.index(not I_p);
         n_n :=nR - n_p;
        P :=zeros(nR, nR);
         for i in 1:nR loop
           P[pivots[i],i] :=1;
         end for;
    //        success :=Modelica.Utilities.Streams.writeRealMatrix(
    //            "U:/OutputX_pTXred_a.mat",
    //            "a",
    //            {a});
          a :=transpose(P)*a;//checked
          ap[1:n_p] :=a[1:n_p];
          an[1:n_n] :=a[n_p + 1:nR];
    //        success :=Modelica.Utilities.Streams.writeRealMatrix(
    //            "U:/OutputX_pTXred_a_.mat",
    //            "a_",
    //            {a});
    //        success :=Modelica.Utilities.Streams.writeRealMatrix(
    //            "U:/OutputX_pTXred_D.mat",
    //            "D",
    //            D);

    //      D :=P*D;
         Dvector :=D*fill(1, nR);
         Dvector :=Dvector[pivots];
         D :=diagonal(Dvector);//checked
         Dp[1:n_p,1:n_p] :=D[1:n_p, 1:n_p];
         Dn[1:n_n,1:n_n] :=D[n_p + 1:nR, n_p + 1:nR];
         Dpvector[1:n_p] :=Dvector[1:n_p];
         Dnvector[1:n_n] :=Dvector[n_p + 1:nR];
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_C.mat",
    //           "C",
    //           C);

         C :=C*P; //checked
         Cp[:,1:n_p] :=C[:, 1:n_p];
         Cn[:,1:n_n] :=C[:, n_p + 1:nR];

    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_C_.mat",
    //           "C_",
    //           C);
    //        success :=Modelica.Utilities.Streams.writeRealMatrix(
    //            "U:/OutputX_pTXred_B.mat",
    //            "B",
    //            B);
         B :=transpose(P)*B; //checked
         Bp[1:n_p,:] :=B[1:n_p, :];
         Bn[1:n_n,:] :=B[n_p + 1:nR, :];
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //            "U:/OutputX_pTXred_B_.mat",
    //            "B_",
    //            B);
         x1[1:nR] :=transpose(P)*x1[1:nR];
         x2[1:nR] :=transpose(P)*x2[1:nR];
    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_x1_.mat",
    //           "x1_",
    //           {x1});

    //      Dvector :=D*fill(1, nR);
          invDvector :=ones(nR) ./ Dvector;
          invDp[1:n_p,1:n_p] :=diagonal(ones(n_p)./Dpvector[1:n_p]);//diagonal(invDvector[1:n_p]);//fill(invDvector[1:n_p], n_p);
          invDn[1:n_n,1:n_n] :=diagonal(ones(n_n)./Dnvector[1:n_n]);//diagonal(invDvector[n_p+1:nR]);//fill(invDvector[n_p + 1:nR], n_n);
    // invD :=fill(invDvector,nR);//diagonal(ones(nR)./Dvector);//{{if i == j then 1/D[i, j] else 0 for i in 1:nR} for j in 1:nR};
    //      invD :=Modelica.Math.Matrices.inv(D);
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_invD.mat",
    //           "invD",
    //           invD);

         J_[1:n_n,1:n_n] :=Dn[1:n_n,1:n_n];//D[n_p + 1:nR, n_p + 1:nR];
         J_[1:n_n,n_n+1:n_n+nX] :=Bn[1:n_n,:];//B[n_p + 1:nR, :];
         J_[n_n+1:n_n+nX,1:n_n] :=Cn[:,1:n_n];//C[:, n_p + 1:nR];
         J_[n_n+1:n_n+nX,n_n+1:n_n+nX] := I - Cp[:,1:n_p]*invDp[1:n_p,1:n_p]*Bp[1:n_p,:];//identity(nX) - C[:, 1:n_p]*(invD[1:n_p, 1:n_p]*B[1:n_p, :]);
    //        J_[nR-n_p+1:n_n+nX,nR-n_p+1:n_n+nX] :=identity(nX) - C[:, 1:n_p]*invDvector[1:n_p]*B[1:n_p, :];

         f_[1:n_n] :=an[1:n_n];//a[n_p + 1:nR];
         f_[n_n+1:n_n+nX] :=b - Cp[:,1:n_p]*invDp[1:n_p,1:n_p]*ap[1:n_p];//b - C[:, 1:n_p]*invD[1:n_p, 1:n_p]*a[1:n_p];
    //       f_[nR-n_p+1:nF-n_p] :=b - C[:, 1:n_p]*(a[1:n_p].*invDvector[1:n_p]);//*(a[1:n_p].*invDvector[1:n_p]);

        // creating Permutation matrix and do pivoting, then solve subproblem with Newton

    //        x2[1:n_p] :=invD[1:n_p, 1:n_p]*(a[1:n_p] - B[1:n_p,:]*x1[nR + 1:nF]);
    //       while not solutionfound1 loop
            x2[n_p+1:nF] :=Initialization.NewtonStep(
             x1[n_p + 1:nF],
             f_[1:n_n+nX],
             J_[1:n_n+nX, 1:n_n+nX],
             n_n+nX);
            x2[1:n_p] :=invDp[1:n_p,1:n_p]*(ap[1:n_p]-Bp[1:n_p,:]*x2[nR+1:nF]);//invD[1:n_p, 1:n_p]*(ap[1:n_p] - B[1:n_p,:]*x2[nR + 1:nF]);

    //          if Modelica.Math.Vectors.norm(abs(x1[n_p+1:nF]-x2[n_p+1:nF]),Modelica.Constants.inf) < eps then
    //             solutionfound1:=true;
    //           end if;
    //          x1[n_p+1:nF] :=x2[n_p + 1:nF];
    //       end while;

    //       for i in 1:n_p loop
    //         x2[i] :=invDvector[i]*(a[i] - B[i, :]*x2[nR + 1:nF]);
    //       end for;
                // x2[1:n_p] :=invD[1:n_p, 1:n_p]*(a[1:n_p] - B[1:n_p,:]*x2[nR + 1:nF]);
    //           x2[1:n_p] :=invDvector[1:n_p].*(a[1:n_p] - B[1:n_p,:]*x2[nR + 1:nF]);

          //   x2[n_p+1:nF] :=Modelica.Math.Matrices.solve(J_[1:nF - n_p, 1:nF - n_p],
          // f_[1:nF - n_p]);

           x2[1:nR] :=P*x2[1:nR];//*transpose(P);
           x1[1:nR] :=P*x1[1:nR];//*transpose(P);

    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //          "U:/OutputX_pTXred_x2.mat",
    //          "x2",
    //          {x2});

    //      x2[1:nR] :=invD*(a - B*x1[nR + 1:nF]);
    //        x2[nR+1:nF] :=Initialization.NewtonStep(
    //        x1[nR + 1:nF],
    //        b - C*invD*a,
    //        I - C*invD*B,
    //        nX);

    //      assert(sum(x2[nR+1:nF]) > 0, "x2[nR+1:nF] is zero at iteration " + String(index));
    //     x2 := P_to_orig_*Initialization.NewtonStep(P_to_id_*x1,f,J,nF);
    //         x2 := Initialization.NewtonStep(x1,f,J,nF);

        x1 :=x2;
         index :=index + 1;

        f :=Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps or index>1000 then
          solutionfound:=true;
        end if;
      end while;
      x2 :=P_to_orig_*x2;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      for i in 1:nF loop
      x_[i] :=max(1e-30,Xfull[i]);
      end for;
    //    x_ :=x2;//Xfull;

       success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_J.mat",
              "J",
              J);

       success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_f.mat",
              "f",
              {f});
      success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_J_.mat",
              "J_",
              J_);

    //  success :=Modelica.Utilities.Streams.writeRealMatrix(
    //         "U:/OutputX_pTXred_x1.mat",
    //         "x1",
    //         {x1});
      //
      // success :=Modelica.Utilities.Streams.writeRealMatrix(
      //        "U:/OutputX_pTXred_x1_.mat",
      //        "x1_",
      //        {x1});

    //   success :=Modelica.Utilities.Streams.writeRealMatrix(
    //     "U:/OutputX_pTXred_B_.mat",
    //     "B_",
    //     B);
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_D_.mat",
         "D_",
         D);
    //
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_n.mat",
         "n",
         {cat(1,pivots,{n_p,n_n})});
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_P.mat",
         "P",
         P);

    end X_pTXred_v4;

    function X_pTXred_v5
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-8;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B_;
      Real[nR] b_;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      MassFraction[nX] Xred_;

      Integer index;
      Boolean success;
      Real[nF,nX] lambda_mass1;

      Real[nR,nR] D;
      Real[nR] Dvector;
      Real[nR] invDvector;
      Real[nR,nX] B;
      Real[nX,nR] C;
      Real[nX,nX] I;
      Real[nR] a;
      Real[nX] b;
      Real[nR,nR] invD;
      Boolean[nR] I_p;
      Integer n_p;
      Integer n_n;
      Integer[nR] pivots;
      Integer[nR,nR] P;
      Real[nF,nF] J_;
      Real[nF] f_;
    algorithm
        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B_[r,firstnonzero[r]] :=1;
           b_[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xred,B_,b_);
         assert(sum(Xfull_in)> 0,"Test");
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;
    //        x1[1:nF] :=cat(1,{0.001*i for i in 1:nF-1},{1-0.001*sum(i for i in 1:nF-1)});//fill(1/nF, nF);
        x1[1:nF] :=fill(1/nF, nF);

         lambda_mass1 :=P_to_orig_*lambda_mass;
         // Xred_ :=Modelica.Math.Matrices.solve(transpose(lambda_mass),Xred);//
         Xred_:=transpose(lambda_mass_)*Xfull_in;

      // Newton algorithm
      solutionfound :=false;
      x1 :=P_to_id_*x1;
       index :=1;
      while not solutionfound loop
    //     f :=Initialization.calc_f.SLE_Tp_v1(x1,Xred,T,p);
    //     J :=Initialization.calc_J.SLE_Tp_v1(x1,T,p);
    //     x2 := Initialization.NewtonStep(x1,f,J,nF);
    //     assert(sum(x1)>0, "x1 is zero at iteration " + String(index));
        f :=-Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps or index>100 then
          solutionfound:=true;
        end if;
         a :=f[1:nR];
         b :=f[nR + 1:nF];
         D :=J[1:nR, 1:nR];
         B :=J[1:nR, nR + 1:nF];
         C :=J[nR + 1:nF, 1:nR];
         I :=identity(nX);
         I_p :={false for i in 1:nR};//{abs(D[i, i]) >= Modelica.Math.Vectors.norm(C[:, i],Modelica.Constants.inf) for i in 1:nR};//

         n_p :=Modelica.Math.BooleanVectors.countTrue(I_p);
         pivots[1:n_p] :=Modelica.Math.BooleanVectors.index(I_p);
         pivots[n_p+1:nR] :=Modelica.Math.BooleanVectors.index(not I_p);
         n_n :=nR - n_p;
        P :=zeros(nR, nR);
         for i in 1:nR loop
           P[pivots[i],i] :=1;
         end for;
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_a.mat",
    //           "a",
    //           {a});
          a :=transpose(P)*a;
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_a_.mat",
    //           "a_",
    //           {a});
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_D.mat",
    //           "D",
    //           D);

    //      D :=P*D;
         Dvector :=D*fill(1, nR);
         Dvector :=Dvector[pivots];
         D :=diagonal(Dvector);
         C :=C*P;
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_B.mat",
    //           "B",
    //           B);
         B :=transpose(P)*B;
    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_x1.mat",
    //           "x1",
    //           {x1});
         x1[1:nR] :=transpose(P)*x1[1:nR];
         x2[1:nR] :=transpose(P)*x2[1:nR];
    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_x1_.mat",
    //           "x1_",
    //           {x1});

    //      Dvector :=D*fill(1, nR);
          invDvector :=ones(nR) ./ Dvector;
         invD :=diagonal(invDvector);//fill(invDvector,nR);//diagonal(ones(nR)./Dvector);//{{if i == j then 1/D[i, j] else 0 for i in 1:nR} for j in 1:nR};
    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //          "U:/OutputX_pTXred_invD.mat",
    //          "invD",
    //          invD);

         J_[1:nR-n_p,1:nR-n_p] :=D[n_p + 1:nR, n_p + 1:nR];
         J_[1:nR-n_p,nR-n_p+1:n_n+nX] :=B[n_p + 1:nR, :];
         J_[nR-n_p+1:n_n+nX,1:nR-n_p] :=C[:, n_p + 1:nR];
            J_[nR-n_p+1:n_n+nX,nR-n_p+1:n_n+nX] :=I - C[:, 1:n_p]*invD[1:n_p, 1:n_p]*B[1:n_p, :]; // CORRECT

               for i in 1:nX loop
                 for j in 1:nX loop
                   J_[nR-n_p+i,nR-n_p+j] :=I[i, j] - C[i, 1:n_p]*(invDvector[1:n_p].*B[1:n_p,j]);
                 end for;
               end for;

    //       J_[nR-n_p+1:n_n+nX,nR-n_p+1:n_n+nX] :=identity(nX) - C[:, 1:n_p]*invDvector[1:n_p]*B[1:n_p, :];

         f_[1:nR-n_p] :=a[n_p + 1:nR];
    //   f_[nR-n_p+1:nF-n_p] :=b - C[:, 1:n_p]*invD[1:n_p, 1:n_p]*a[1:n_p]; // CORRECT
               for i in 1:nX loop
                  f_[nR-n_p+i] :=b[i] - C[i, 1:n_p]*(invDvector[1:n_p].*a[1:n_p]);//invD[1:n_p,1:n_p]*a[1:n_p];
    //              f_[nR-n_p+i] :=b[i] - C[i, 1:n_p]*invD[1:n_p,1:n_p]*a[1:n_p];
               end for;
    //       f_[nR-n_p+1:nF-n_p] :=b - C[:, 1:n_p]*(a[1:n_p].*invDvector[1:n_p]);//*(a[1:n_p].*invDvector[1:n_p]);

        // creating Permutation matrix and do pivoting, then solve subproblem with Newton

    //        x2[1:n_p] :=invD[1:n_p, 1:n_p]*(a[1:n_p] - B[1:n_p]*x1[nR + 1:nF]);
    //        while Modelica.Math.Vectors.norm(x2[n_p+1:nF]-x1[n_p+1:nF],Modelica.Constants.inf) > eps loop
           x2[n_p+1:nF] :=Initialization.NewtonStep(
            x1[n_p + 1:nF],
            f_[1:nF - n_p],
            J_[1:nF - n_p, 1:nF - n_p],
            nF - n_p);
    //          x1[n_p+1:nF] :=x2[n_p + 1:nF];
    //        end while;
    //         for i in 1:n_p loop
    //           x2[i] :=invD[i, i]*(a[i] - B[i, :]*x2[nR + 1:nF]);
    //         end for;

              x2[1:n_p] :=invD[1:n_p, 1:n_p]*(a[1:n_p] - B[1:n_p,:]*x2[nR + 1:nF]);//CORRECT

    //         x2[1:n_p] :=invDvector[1:n_p].*(a[1:n_p] - B[1:n_p]*x2[nR + 1:nF]);

           x2[1:nR] :=P*x2[1:nR];//*transpose(P);
    //        x1[1:nR] :=P*x1[1:nR];//*transpose(P);

    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //          "U:/OutputX_pTXred_x2.mat",
    //          "x2",
    //          {x2});

    //      x2[1:nR] :=invD*(a - B*x1[nR + 1:nF]);
    //        x2[nR+1:nF] :=Initialization.NewtonStep(
    //        x1[nR + 1:nF],
    //        b - C*invD*a,
    //        I - C*invD*B,
    //        nX);

    //      assert(sum(x2[nR+1:nF]) > 0, "x2[nR+1:nF] is zero at iteration " + String(index));
    //     x2 := P_to_orig_*Initialization.NewtonStep(P_to_id_*x1,f,J,nF);
    //         x2 := Initialization.NewtonStep(x1,f,J,nF);

        x1 :=x2;
         index :=index + 1;
      end while;
      x2 :=P_to_orig_*x2;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      for i in 1:nF loop
      x_[i] :=max(1e-20,Xfull[i]);
      end for;
       x_ :=Xfull;

       success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_J.mat",
              "J",
              J);

       success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_f.mat",
              "f",
              {f});

       success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_f_.mat",
              "f_",
              {f_});
      success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_J_.mat",
              "J_",
              J_);

    //  success :=Modelica.Utilities.Streams.writeRealMatrix(
    //         "U:/OutputX_pTXred_x1.mat",
    //         "x1",
    //         {x1});
      //
      // success :=Modelica.Utilities.Streams.writeRealMatrix(
      //        "U:/OutputX_pTXred_x1_.mat",
      //        "x1_",
      //        {x1});

       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_invD.mat",
         "invD",
         invD);
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_D_.mat",
         "D_",
         D);
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_Xred.mat",
         "Xred",
         {Xred_});
    //    success :=Modelica.Utilities.Streams.writeRealMatrix(
    //      "U:/OutputX_pTXred_b.mat",
    //      "b",
    //      {b});
    //
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_n.mat",
         "n",
         {cat(1,pivots,{n_p,n_n},b)});

       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_P_nu.mat",
         "P",
         P_to_id_);

       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_lambda_id.mat",
         "lambda_id",
         lambda_id);

       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_lambda_.mat",
         "lambda_",
         lambda_);
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_P.mat",
         "P",
         P);

    end X_pTXred_v5;
  end MixtureLiquid;

  package MixtureLiquidAlgorithm
    import Modelica.Math;
    import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;

    extends Modelica.Media.Interfaces.PartialMixtureMedium(
       ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX,
       reference_X=userInterface.refX,
       nX = nS-nR,
       reducedX = true,
       nXi = nX-1,
       mediumName = userInterface.mediumName,
       substanceNames=cat(1,userInterface.solidSubstanceNames,userInterface.liquidSubstanceNames),
       Density(start=10, nominal=10),
       AbsolutePressure(start=10e5, nominal=10e5),
       Temperature(min=273.15, max=573.15, start=300, nominal=300),
       T_default = userInterface.Tstart,
       p_default = userInterface.pstart);
  //      nXi = nX-1,
  //      fixedX=false,
  //      reducedX = true,
  //      singleState=false,

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

      Real[nR,ns+nL] nu = zeros(nR,ns+nL) "Stoichiometry matrix of solid-liquid and dissociation equilibria" annotation(Dialog(tab="Reaction"));

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
    constant Real[nR,nF] nu_id_= Media.Common.Reaction.calc_nu_id(nu);
                                                                      //nu
    constant Integer[nF,nF] P_to_orig_=Media.Common.Reaction.calc_P_nu_id(nu);
                                                                              //identity(nF);
    constant Integer[nF,nF] P_to_id_ = transpose(P_to_orig_);
    constant Real[nR,nF] nu_id__ = nu_id_*transpose(P_to_orig_);
    constant Real[nF,nX] lambda_nu = transpose(cat(2,-transpose(nu_id_[:,nR+1:nF]),identity(nX)));
  //   constant Real[nF,nX] lambda_id = cat(2,identity(nX),transpose(nu_id[1:nX,nX+1:nF]));
    constant Real[nR,nF] nu_mass = Functions.Reaction.calc_nu_mass(nu);
    constant Real[nF] MMX_ = P_to_id_*MMX;

    constant UserInterface
      userInterface "User interface"
      annotation (Placement(transformation(extent={{-10,-8},{10,12}})));

    constant Real[nF,nX] lambda_=
      Media.Common.Reaction.calc_lambda(
      nu,
      nR,
      nF);
    constant Real[nF,nX] lambda_id = Media.Common.Reaction.calc_lambda_id(lambda_);
    constant Real[nF,nX] lambda = P_to_orig*lambda_id;
    constant Integer[nF,nF] P_to_orig= Media.Common.Reaction.calc_P_lambda_id(lambda_);
    constant Integer[nF,nF] P_to_id = transpose(P_to_orig);
    constant Real[nR,nF] nu_id = nu*transpose(P_to_id);
    constant Real[nF,nX] lambda_mass_=P_to_orig_*{{lambda_nu[i,j]/MMX_[i] for j in 1:nX} for i in 1:nF} "Null space of mass based stoichiometry matrix";
    constant Real[nF,nX] lambda_mass={{lambda_[i,j]/MMX[i] for j in 1:nX} for i in 1:nF} "Null space of mass based stoichiometry matrix";
  //  constant Real[nX] scale = {sum(lambda_mass_[i,j] for i in 1:nF) for j in 1:nX};
  //   constant Real[:,:] lambda_mass = {{lambda_mass_[i,j]/scale[j] for j in 1:nX} for i in 1:nF};

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
      //final standardOrderComponents=true)
      MoleFraction[ns+nL] Yfull(start=Yfullstart) "Full mole fraction vector";
      MassFraction[ns+nL] Xfull(start=Xfullstart,each min = 0) "Full mass fraction vector";
    //   SI.Mass[ns+nL] mfull(start=mfullstart*scale)  "Full amount of substance vector (substitute)";
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
      //     Real[ns+nL] logreacBase(start = cat(1,-zstart[1:ns],logastart));//-zstart[ns+1:ns+nL]));
      Real pH;
      Real I;
      //     Real z [ns](start = zstart[1:ns]);
      MoleFraction[ns] Ys;
      MoleFraction[nL] Yl;
      MassFraction[ns] Xs;
      MassFraction[nL] Xl;
      SI.VolumeFraction Phil;

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
      parameter Real scale = 100000;

      //Real[nX] X_ "Substitute for X to ensure non-reactive species to be > 0";

      //   initial algorithm
      //      Xfull :=Xfullstart;

    algorithm
      //     when Modelica.Math.Vectors.norm(X-pre(X),Modelica.Constants.inf) > 1e-4 then
      Xfull :=X_pTXred(
        p,
        T,
        X,Xfull);
      //     end when;

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

    //   mred_ = transpose(lambda_mass)*mfull;
    //   mred_ = {max(1e-25,mred[i]) for i in 1:nX};
    //   X = mred_*scale;

    //   Xfull[1:nF] = mfull[1:nF]/sum(mfull);
    //   X_ = transpose(lambda_mass_scaled)*mfull;
    //   X_ = {max(1e-25,X[i]) for i in 1:nX};

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
      //     logreacBase = cat(1,-z[1:ns],loga);//-z[1+ns:ns+nL]);
      //     nu*logreacBase = logK;
      a = exp(loga);
      pH = -loga[Hindex-ns]*log10(exp(1));
      Phil = Functions.calc_Phil(T,p,Xfull);

      //Fischer-Burmeister function
    //    ones(ns)*eps = Xfull[1:ns].*z;
      // ones(ns+nL)*1e-20 = mfull.*z;
    //    zeros(ns+nL) = Xfull + z - sqrt(Xfull.^2+z.^2 + ones(ns+nL)*1e-20);
      //       zeros(ns) = Xfull[1:ns] + z - sqrt(Xfull[1:ns].^2+z.^2 .+ eps);

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
      Real[nF + ns] init=X_pTXred(
            p,
            T,
            Xred);
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
      parameter Real init = 1e-10;
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
      Real[nF,nX] lambda_mass_id = P_to_id_ * lambda_mass;
    algorithm
      //ensure no species amount to be zero for numerical issues
      for i in 1:nX loop
        Xrednonzero[i] :=max(Modelica.Constants.eps, Xred[i]);
      end for;

      Xs :=Modelica.Math.Matrices.solve(transpose(lambda_mass_id[nR+1:nF, :]), Xrednonzero);
      Xfull :=cat(1,zeros(nF - nX),Xs);
      Xfull :=P_to_orig_*Xfull;

      Xred_ :=transpose(lambda_mass_)*Xfull;

      // find initial guess x1 for Newton solver
      if sum(Xinit)> 0 then
      // ninit from initial guess provided as input (Xinit)
         scale :=sum(transpose(lambda_mass_)*Xinit);
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
      x1 :=P_to_id_*x1;
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
      x2 :=P_to_orig_*x2;
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

      Real[nF+ns] X_;

    algorithm
      //calculate start values for inner-outer algorithm
      T :=exp((log(Tmin) + log(Tmax))/2);
      X_ :=X_pTXred(
          p,
          T,
          Xred);
      X :=X_[1:nF];

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
        X_ :=X_pTXred(
            p,
            T,
            Xred,
            X_);
        X :=X_[1:nF];
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

      Real[nF+ns] X_;

    algorithm
      //calculate start values for inner-outer algorithm
      T :=exp((log(Tmin) + log(Tmax))/2);
      X_ :=X_pTXred(
          p,
          T,
          Xred);
      X :=X_[1:nF];

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
        X_ :=X_pTXred(
            p,
            T,
            Xred,
            X_);
        X :=X_[1:nF];
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

      Real[nF+ns] X_;
    algorithm
      //calculate start values for inner-outer algorithm
      p :=exp((log(pmin) + log(pmax))/2);//1.1*pmin;//
      X_ :=X_pTXred(
          p,
          T,
          Xred);
      X :=X_[1:nF];

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
        X_ :=X_pTXred(
            p,
            T,
            Xred,
            X_);
        X :=X_[1:nF];
      end while;
    end Xp_dTXred;

    function calc_Xred
      input SI.MassFraction[ns+nL] Xfull;
      output SI.MassFraction[nX] Xred;

    algorithm
      Xred :=transpose(lambda_mass)*Xfull/(sum(transpose(lambda_mass)*Xfull));

    end calc_Xred;

    package Functions
      extends Common.Functions(nSfunfun= ns, dataSfunfun = datas, nLfunfun = nL, dataLfunfun = datal, LiquidModelfunfun = LiquidModel, interactionLfunfun = interactionL);
    end Functions;

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

          x_orig :=P_to_orig_*x;
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
          gi_id :=P_to_id_*gi;

        //    gr :=nu_id*gi_id;
           gr :=nu_id__*gi;
          logK :=-gr/Modelica.Constants.R/T;

          a :=gamma_i .* m;
          logreacBase :=cat(1,-z[1:ns],log(a)-z[ns+1:ns+nL]);
           logreacBase_id :=P_to_id_*logreacBase;

        //   x_id :=P_to_id_*x;

          //isopotential
           f[1:nR] :=nu_id__*logreacBase - logK;
        //   f[1:nR] :=nu_id*logreacBase_id - logK;

          //reduced mass balance
          // f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
        //   f[nR+1:nF] :=transpose(lambda_nu)*x_id - Xred;//transpose(lambda_mass)*mass - Xred;
          f[nR+1:nF] :=transpose(lambda_nu)*x - Xred;//transpose(lambda_mass)*mass - Xred;
        end SLE_Tp;

        function SLE_Tp_
          "SLE initialization according to Leal et al. (2016): calculation of residuals"
          input Real[nF] x;
          input Real[nX] Xred;
          input SI.Temperature T;
          input SI.Pressure p;
          output Real[nF] f;

        protected
          SI.MolarEnergy gi[ns+nL];
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
          Real[ns] slack = 1e-20./x[1:ns];//x[nF+1:nF+ns];
          parameter Real eps = 1e-15;

          Real[ns+nL] logreacBase;

          SI.Mass[ns+nL] mass;

        algorithm
          //assert(sum(x[1:ns+nL]) > 0, "x is zero in SLE_Tp");

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          Xs :=Functions.calc_X_M(Ys, MMX[1:ns]);
          Xl :=Functions.calc_X_M(Yl, MMX[1 + ns:ns + nL]);

          m[1:nL-1] := Functions.LiquidFunctions.calc_mfromY(Yl);
          m[nL] := Yl[nL]/MH2O;

          gamma_i := Functions.LiquidFunctions.calc_gamma(T,p,Xl);
          gi := Functions.calc_gi(T, p) .* MMX;

          gr :=nu*gi;
          logK :=-gr/Modelica.Constants.R/T;

          a :=gamma_i .* m;
          logreacBase :=cat(1,-slack,log(a));

          mass :=x[1:nF] .* MMX;

          //isopotential
          f[1:nR] :=nu*logreacBase - logK;

          //reduced mass balance
          f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
        end SLE_Tp_;

        function SLE_Tp_v1
          "SLE initialization according to Leal et al. (2016): calculation of residuals"
          input Real[nF] x;
          input Real[nX] Xred;
          input SI.Temperature T;
          input SI.Pressure p;
          output Real[nF] f;

        protected
          SI.MolarEnergy gi[ns+nL];
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
          Real tau = 1e-30;
          SI.MassFraction[ns+nL] z = {if x[i] > 0 then tau/max(tau,x[i]) else 0 for i in 1:ns+nL};//tau./x;

          Real[ns+nL] logreacBase;

          SI.Mass[ns+nL] mass;

        algorithm
          assert(sum(x[1:ns+nL]) > 0, "x is zero in SLE_Tp");

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          Xs :=Functions.calc_X_M(Ys, MMX[1:ns]);
          Xl :=Functions.calc_X_M(Yl, MMX[1 + ns:ns + nL]);

          m[1:nL-1] := Functions.LiquidFunctions.calc_mfromY(Yl);
          m[nL] := Yl[nL]/MH2O;

          gamma_i := Functions.LiquidFunctions.calc_gamma(T,p,Xl);
          gi := Functions.calc_gi(T, p) .* MMX;

          gr :=nu*gi;
          logK :=-gr/Modelica.Constants.R/T;

          a :=gamma_i .* m;
          logreacBase :=cat(1,-z[1:ns],log(a)-z[ns+1:ns+nL]);

          mass :=x[1:nF] .* MMX;

          //isopotential
          f[1:nR] :=nu*logreacBase - logK;

          //reduced mass balance
          // f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
          f[nR+1:nF] :=transpose(lambda_)*x - Xred;//transpose(lambda_mass)*mass - Xred;
        end SLE_Tp_v1;

        function SLE_Tp_v2
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
          Real tau = 1e-30;
          SI.MassFraction[ns+nL] z = {if x[i] > 0 then tau/max(tau,x[i]) else 0 for i in 1:ns+nL};//tau./x;

          Real[ns+nL] logreacBase;
          Real[ns+nL] logreacBase_id;
          Real[ns+nL] x_id;

          SI.Mass[ns+nL] mass;

        algorithm
          assert(sum(x[1:ns+nL]) > 0, "x is zero in SLE_Tp");

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
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

          mass :=x[1:nF] .* MMX;
          x_id :=P_to_id*x;

          //isopotential
           f[1:nR] :=nu*logreacBase - logK;
        //   f[1:nR] :=nu_id*logreacBase_id - logK;

          //reduced mass balance
          // f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
          f[nR+1:nF] :=transpose(lambda_id)*x_id - Xred;//transpose(lambda_mass)*mass - Xred;
        end SLE_Tp_v2;

        function SLE_Tp_v3
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
          Real tau = 1e-30;
          SI.MassFraction[ns+nL] z;// = {if x[i] > 0 then tau/max(tau,x[i]) else 0 for i in 1:ns+nL};//tau./x;

          Real[ns+nL] logreacBase;
          Real[ns+nL] logreacBase_id;
          Real[ns+nL] x_orig;

          SI.Mass[ns+nL] mass;

        algorithm

          assert(sum(x[1:ns+nL]) > 0, "x is zero in SLE_Tp");

          x_orig :=P_to_orig_*x;

          z :={if x_orig[i] > 0 then tau/max(tau, x_orig[i]) else 0 for i in 1:ns + nL};

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
          gi_id :=P_to_id_*gi;

        //    gr :=nu_id*gi_id;
           gr :=nu_id__*gi;
          logK :=-gr/Modelica.Constants.R/T;

          a :=gamma_i .* m;
          logreacBase :=cat(1,-z[1:ns],log(a)-z[ns+1:ns+nL]);
           logreacBase_id :=P_to_id_*logreacBase;

        //   x_id :=P_to_id_*x;

          //isopotential
           f[1:nR] :=nu_id__*logreacBase - logK;
        //   f[1:nR] :=nu_id*logreacBase_id - logK;

          //reduced mass balance
          // f[nR+1:nF] :=transpose(lambda_mass)*mass - Xred;
        //   f[nR+1:nF] :=transpose(lambda_nu)*x_id - Xred;//transpose(lambda_mass)*mass - Xred;
          f[nR+1:nF] :=transpose(lambda_nu)*x - Xred;//transpose(lambda_mass)*mass - Xred;
        end SLE_Tp_v3;
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

          x_orig :=P_to_orig_*x;

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
          H_id :=transpose(P_to_id_)*H;

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
             J[1:nR,:] :=nu_id__*H;//nu_id__ .* H;//nu_id__*H;//
             J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id_);//transpose(P_nu_id);

        //      J[1:nR,:] :=nu_id_ .* H_id;
        //     J[1:nR,:] :=nu_id .* H_id;

          //reduced mass balance
           J[nR+1:nF,:] :=transpose(lambda_nu);
        end SLE_Tp;

        function SLE_Tp_
          "SLE initialization according to Leal et al. (2016): calculation of Jacobian"
          input Real[nF] x;
          input SI.Temperature T;
          input SI.Pressure p;
          output Real[nF,nF] J;

        protected
          SI.MolarEnergy gi[ns+nL];
          SI.MoleFraction[ns+nL] Y;
          SI.MassFraction[ns+nL] X;
          SI.MoleFraction[nL] Yl;
          SI.MassFraction[nL] Xl;
          SI.MoleFraction[ns] Ys;
          SI.MassFraction[ns] Xs;
          parameter Real eps = 1e-15;
          Real[ns] slack = 1e-20./x[1:ns];//x[nF+1:nF+ns];
          Real[ns+nL] logreacBase;
        algorithm

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          for i in 1:ns+nL loop
            //isopotential
            for r in 1:nR loop
              if i < ns+1 then
                J[r,i] := if x[i] > 0 then nu[r,i]*slack[i]/x[i] else 0;
               else
                 J[r,i] :=if nu[r, i] > 0 then nu[r, i]*(1 - Yl[i - ns])/x[i] else 0;
               end if;
            end for;

            //reduced mass balance
            for j in 1:ns+nL - nR loop
              J[j + nR,i] := lambda_mass[i, j]*MMX[i];
            end for;
          end for;
        end SLE_Tp_;

        function SLE_Tp_v1
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
          Real tau = 1e-20;
          SI.MassFraction[ns+nL] z;
          Real[ns+nL,ns+nL] H;
          Real[ns+nL] diagH;
        algorithm

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          for i in 1:ns+nL loop
            n[i] :=max(tau, x[i]);
          end for;
          z :=tau./n;

          diagH[1:ns] :=tau ./ n[1:ns] .^ 2;
          diagH[ns+1:ns+nL] :=(1 .- Yl) ./ n[ns + 1:ns + nL] + tau ./ n[ns+1:ns+nL] .^ 2;

          H :=fill(diagH,ns+nL);

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
          J[1:nR,:] :=nu .* H;

          //reduced mass balance
           J[nR+1:nF,:] :=transpose(lambda_);
        end SLE_Tp_v1;

        function SLE_Tp_v2
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
          Real tau = 1e-20;
          SI.MassFraction[ns+nL] z;
          Real[ns+nL,ns+nL] H;
          Real[ns+nL,ns+nL] H_id;
          Real[ns+nL] diagH;
        algorithm

          Y :=x[1:ns+nL]/sum(x[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          for i in 1:ns+nL loop
            n[i] :=max(tau, x[i]);
          end for;
          z :=tau./n;

          diagH[1:ns] :=tau ./ n[1:ns] .^ 2;
          diagH[ns+1:ns+nL] :=(1 .- Yl) ./ n[ns + 1:ns + nL] + tau ./ n[ns+1:ns+nL] .^ 2;

          H :=fill(diagH,ns+nL);
          H_id :=P_to_id*H;

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
             J[1:nR,:] :=nu .* H;
             J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id);
        //     J[1:nR,:] :=nu_id .* H_id;

          //reduced mass balance
           J[nR+1:nF,:] :=transpose(lambda_id);
        end SLE_Tp_v2;

        function SLE_Tp_v3
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
          Real tau = 1e-20;
          SI.MassFraction[ns+nL] z;
          Real[ns+nL,ns+nL] H;
          Real[ns+nL,ns+nL] H_id;
          Real[ns+nL] diagH;
          Real[ns+nL] x_orig;
        algorithm

          x_orig :=P_to_orig_*x;

          Y :=x_orig[1:ns+nL]/sum(x_orig[1:ns+nL]);
          X :=Functions.calc_Xfull(Y);

          Ys :=Y[1:ns]/sum(Y[1:ns]);
          Yl :=Y[1 + ns:ns + nL]/sum(Y[1 + ns:ns + nL]);

          for i in 1:ns+nL loop
            n[i] :=max(tau, x_orig[i]);
          end for;
          z :=tau./n;

          diagH[1:ns] :=tau ./ n[1:ns] .^ 2;
          diagH[ns+1:ns+nL] :=(1 .- Yl) ./ n[ns + 1:ns + nL] + tau ./ n[ns+1:ns+nL] .^ 2;

          H :=fill(diagH,ns+nL);
          H_id :=P_to_id_*H;

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
             J[1:nR,:] :=nu_id__ .* H;
             J[1:nR,:] :=J[1:nR, :]*transpose(P_to_id_);//transpose(P_nu_id);
        //     J[1:nR,:] :=nu_id .* H_id;

          //reduced mass balance
           J[nR+1:nF,:] :=transpose(lambda_nu);
        end SLE_Tp_v3;
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

    function X_pTXred_
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-8;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B;
      Real[nR] b;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      Integer index;
    algorithm
        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B[r,firstnonzero[r]] :=1;
           b[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xred,B,b);
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;

      // Newton algorithm
      solutionfound :=false;
      while not solutionfound loop
        f :=Initialization.calc_f.SLE_Tp(x1,Xred,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);
        x2 := Initialization.NewtonStep(x1,f,J,nF);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
          solutionfound:=true;
        end if;
        x1 :=x2;
      end while;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      for i in 1:nF loop
      x_[i] :=max(1e-20,Xfull[i]);
      end for;

    end X_pTXred_;

    function X_pTXred_v2
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-8;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B;
      Real[nR] b;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      MassFraction[nX] Xred_;

      Integer index;
    algorithm
        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B[r,firstnonzero[r]] :=1;
           b[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xred,B,b);
         assert(sum(Xfull_in)> 0,"Test");
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;

         Xred_ :=transpose(lambda_mass_)*Xfull_in;
      // Newton algorithm
      solutionfound :=false;
      while not solutionfound loop
    //     f :=Initialization.calc_f.SLE_Tp_v1(x1,Xred,T,p);
    //     J :=Initialization.calc_J.SLE_Tp_v1(x1,T,p);
    //     x2 := Initialization.NewtonStep(x1,f,J,nF);

        f :=Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);
        x2 := P_to_orig_*Initialization.NewtonStep(P_to_id_*x1,f,J,nF);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
          solutionfound:=true;
        end if;
        x1 :=x2;
      end while;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      for i in 1:nF loop
      x_[i] :=max(1e-20,Xfull[i]);
      end for;

    end X_pTXred_v2;

    function X_pTXred_v3
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-8;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B_;
      Real[nR] b_;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      MassFraction[nX] Xred_;

      Integer index;
      Boolean success;
      Real[nF,nX] lambda_mass1;

      Real[nR,nR] D;
      Real[nR,nX] B;
      Real[nX,nR] C;
      Real[nX,nX] I;
      Real[nR] a;
      Real[nX] b;
      Real[nR,nR] invD;
      Integer nNS;
    algorithm
        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B_[r,firstnonzero[r]] :=1;
           b_[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xred,B_,b_);
         assert(sum(Xfull_in)> 0,"Test");
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;

         lambda_mass1 :=P_to_orig_*lambda_mass;
         // Xred_ :=Modelica.Math.Matrices.solve(transpose(lambda_mass),Xred);//
         Xred_:=transpose(lambda_mass_)*Xfull_in;
         success :=Modelica.Utilities.Streams.writeRealMatrix(
        "U:/OutputX_pTXred.mat",
        "J",
        lambda_mass1);
      // Newton algorithm
      solutionfound :=false;
      x1 :=P_to_id_*x1;
      while not solutionfound loop
    //     f :=Initialization.calc_f.SLE_Tp_v1(x1,Xred,T,p);
    //     J :=Initialization.calc_J.SLE_Tp_v1(x1,T,p);
    //     x2 := Initialization.NewtonStep(x1,f,J,nF);

        f :=Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);
        a :=f[1:nR];
        b :=f[nR + 1:nF];
        D :=J[1:nR, 1:nR];
        B :=J[1:nR, nR + 1:nF];
        C :=J[nR + 1:nF, 1:nR];
        I :=identity(nX);
        invD :={{if i == j then 1/D[i, j] else 0 for i in 1:nR} for j in 1:nR};
        x2[nR+1:nF] :=Initialization.NewtonStep(
          x1[nR + 1:nF],
          b - C*invD*a,
          I - C*invD*B,
          nX);
        x2[1:nR] :=invD*(a - B*x2[nR + 1:nF]);
    //     x2 := P_to_orig_*Initialization.NewtonStep(P_to_id_*x1,f,J,nF);
         x2 := Initialization.NewtonStep(x1,f,J,nF);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps then
          solutionfound:=true;
        end if;
        x1 :=x2;
      end while;
      x2 :=P_to_orig_*x2;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      for i in 1:nF loop
      x_[i] :=max(1e-20,Xfull[i]);
      end for;

    end X_pTXred_v3;

    function X_pTXred_v4
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-3;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;
      Boolean solutionfound1 = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B_;
      Real[nR] b_;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      MassFraction[nX] Xred_;

      Integer index;
      Boolean success;
      Real[nF,nX] lambda_mass1;

      Real[nR,nR] D;
      Real[nR,nR] Dp;
      Real[nR,nR] Dn;
      Real[nR,nR] invDp;
      Real[nR,nR] invDn;
      Real[nR] Dvector;
      Real[nR] Dpvector;
      Real[nR] Dnvector;
      Real[nR] invDvector;
      Real[nR,nX] B;
      Real[nR,nX] Bp;
      Real[nR,nX] Bn;
      Real[nX,nR] C;
      Real[nX,nX] I;
      Real[nR] a;
      Real[nX] b;
      Real[nR,nR] invD;
      Boolean[nR] I_p;
      Integer n_p;
      Integer n_n;
      Integer[nR] pivots;
      Integer[nR,nR] P;
      Real[nF,nF] J_;
      Real[nF] f_;

      Real[nR] ap;
      Real[nR] an;
      Real[nX,nR] Cp;
      Real[nX,nR] Cn;
    algorithm
        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B_[r,firstnonzero[r]] :=1;
           b_[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xred,B_,b_);
         assert(sum(Xfull_in)> 0,"Test");
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;
    //        x1[1:nF] :=cat(1,{0.001*i for i in 1:nF-1},{1-0.001*sum(i for i in 1:nF-1)});//fill(1/nF, nF);
        x1[1:nF] :=fill(1/nF, nF);

         lambda_mass1 :=P_to_orig_*lambda_mass;
         // Xred_ :=Modelica.Math.Matrices.solve(transpose(lambda_mass),Xred);//
         Xred_:=transpose(lambda_mass_)*Xfull_in;

      // Newton algorithm
      solutionfound :=false;
      x1 :=P_to_id_*x1;
       index :=1;
        f :=Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);
      while not solutionfound loop
    //     f :=Initialization.calc_f.SLE_Tp_v1(x1,Xred,T,p);
    //     J :=Initialization.calc_J.SLE_Tp_v1(x1,T,p);
    //     x2 := Initialization.NewtonStep(x1,f,J,nF);
    //     assert(sum(x1)>0, "x1 is zero at iteration " + String(index));

         a :=f[1:nR];
         b :=f[nR + 1:nF];
         D :=J[1:nR, 1:nR];
         B :=J[1:nR, nR + 1:nF];
         C :=J[nR + 1:nF, 1:nR];
         I :=identity(nX);
         I_p :={false for i in 1:nR};//{abs(D[i, i]) >= Modelica.Math.Vectors.norm(C[:, i],Modelica.Constants.inf) for i in 1:nR};

         n_p :=Modelica.Math.BooleanVectors.countTrue(I_p);
         pivots[1:n_p] :=Modelica.Math.BooleanVectors.index(I_p);
         pivots[n_p+1:nR] :=Modelica.Math.BooleanVectors.index(not I_p);
         n_n :=nR - n_p;
        P :=zeros(nR, nR);
         for i in 1:nR loop
           P[pivots[i],i] :=1;
         end for;
    //        success :=Modelica.Utilities.Streams.writeRealMatrix(
    //            "U:/OutputX_pTXred_a.mat",
    //            "a",
    //            {a});
          a :=transpose(P)*a;//checked
          ap[1:n_p] :=a[1:n_p];
          an[1:n_n] :=a[n_p + 1:nR];
    //        success :=Modelica.Utilities.Streams.writeRealMatrix(
    //            "U:/OutputX_pTXred_a_.mat",
    //            "a_",
    //            {a});
    //        success :=Modelica.Utilities.Streams.writeRealMatrix(
    //            "U:/OutputX_pTXred_D.mat",
    //            "D",
    //            D);

    //      D :=P*D;
         Dvector :=D*fill(1, nR);
         Dvector :=Dvector[pivots];
         D :=diagonal(Dvector);//checked
         Dp[1:n_p,1:n_p] :=D[1:n_p, 1:n_p];
         Dn[1:n_n,1:n_n] :=D[n_p + 1:nR, n_p + 1:nR];
         Dpvector[1:n_p] :=Dvector[1:n_p];
         Dnvector[1:n_n] :=Dvector[n_p + 1:nR];
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_C.mat",
    //           "C",
    //           C);

         C :=C*P; //checked
         Cp[:,1:n_p] :=C[:, 1:n_p];
         Cn[:,1:n_n] :=C[:, n_p + 1:nR];

    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_C_.mat",
    //           "C_",
    //           C);
    //        success :=Modelica.Utilities.Streams.writeRealMatrix(
    //            "U:/OutputX_pTXred_B.mat",
    //            "B",
    //            B);
         B :=transpose(P)*B; //checked
         Bp[1:n_p,:] :=B[1:n_p, :];
         Bn[1:n_n,:] :=B[n_p + 1:nR, :];
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //            "U:/OutputX_pTXred_B_.mat",
    //            "B_",
    //            B);
         x1[1:nR] :=transpose(P)*x1[1:nR];
         x2[1:nR] :=transpose(P)*x2[1:nR];
    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_x1_.mat",
    //           "x1_",
    //           {x1});

    //      Dvector :=D*fill(1, nR);
          invDvector :=ones(nR) ./ Dvector;
          invDp[1:n_p,1:n_p] :=diagonal(ones(n_p)./Dpvector[1:n_p]);//diagonal(invDvector[1:n_p]);//fill(invDvector[1:n_p], n_p);
          invDn[1:n_n,1:n_n] :=diagonal(ones(n_n)./Dnvector[1:n_n]);//diagonal(invDvector[n_p+1:nR]);//fill(invDvector[n_p + 1:nR], n_n);
    // invD :=fill(invDvector,nR);//diagonal(ones(nR)./Dvector);//{{if i == j then 1/D[i, j] else 0 for i in 1:nR} for j in 1:nR};
    //      invD :=Modelica.Math.Matrices.inv(D);
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_invD.mat",
    //           "invD",
    //           invD);

         J_[1:n_n,1:n_n] :=Dn[1:n_n,1:n_n];//D[n_p + 1:nR, n_p + 1:nR];
         J_[1:n_n,n_n+1:n_n+nX] :=Bn[1:n_n,:];//B[n_p + 1:nR, :];
         J_[n_n+1:n_n+nX,1:n_n] :=Cn[:,1:n_n];//C[:, n_p + 1:nR];
         J_[n_n+1:n_n+nX,n_n+1:n_n+nX] := I - Cp[:,1:n_p]*invDp[1:n_p,1:n_p]*Bp[1:n_p,:];//identity(nX) - C[:, 1:n_p]*(invD[1:n_p, 1:n_p]*B[1:n_p, :]);
    //        J_[nR-n_p+1:n_n+nX,nR-n_p+1:n_n+nX] :=identity(nX) - C[:, 1:n_p]*invDvector[1:n_p]*B[1:n_p, :];

         f_[1:n_n] :=an[1:n_n];//a[n_p + 1:nR];
         f_[n_n+1:n_n+nX] :=b - Cp[:,1:n_p]*invDp[1:n_p,1:n_p]*ap[1:n_p];//b - C[:, 1:n_p]*invD[1:n_p, 1:n_p]*a[1:n_p];
    //       f_[nR-n_p+1:nF-n_p] :=b - C[:, 1:n_p]*(a[1:n_p].*invDvector[1:n_p]);//*(a[1:n_p].*invDvector[1:n_p]);

        // creating Permutation matrix and do pivoting, then solve subproblem with Newton

    //        x2[1:n_p] :=invD[1:n_p, 1:n_p]*(a[1:n_p] - B[1:n_p,:]*x1[nR + 1:nF]);
    //       while not solutionfound1 loop
            x2[n_p+1:nF] :=Initialization.NewtonStep(
             x1[n_p + 1:nF],
             f_[1:n_n+nX],
             J_[1:n_n+nX, 1:n_n+nX],
             n_n+nX);
            x2[1:n_p] :=invDp[1:n_p,1:n_p]*(ap[1:n_p]-Bp[1:n_p,:]*x2[nR+1:nF]);//invD[1:n_p, 1:n_p]*(ap[1:n_p] - B[1:n_p,:]*x2[nR + 1:nF]);

    //          if Modelica.Math.Vectors.norm(abs(x1[n_p+1:nF]-x2[n_p+1:nF]),Modelica.Constants.inf) < eps then
    //             solutionfound1:=true;
    //           end if;
    //          x1[n_p+1:nF] :=x2[n_p + 1:nF];
    //       end while;

    //       for i in 1:n_p loop
    //         x2[i] :=invDvector[i]*(a[i] - B[i, :]*x2[nR + 1:nF]);
    //       end for;
                // x2[1:n_p] :=invD[1:n_p, 1:n_p]*(a[1:n_p] - B[1:n_p,:]*x2[nR + 1:nF]);
    //           x2[1:n_p] :=invDvector[1:n_p].*(a[1:n_p] - B[1:n_p,:]*x2[nR + 1:nF]);

          //   x2[n_p+1:nF] :=Modelica.Math.Matrices.solve(J_[1:nF - n_p, 1:nF - n_p],
          // f_[1:nF - n_p]);

           x2[1:nR] :=P*x2[1:nR];//*transpose(P);
           x1[1:nR] :=P*x1[1:nR];//*transpose(P);

    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //          "U:/OutputX_pTXred_x2.mat",
    //          "x2",
    //          {x2});

    //      x2[1:nR] :=invD*(a - B*x1[nR + 1:nF]);
    //        x2[nR+1:nF] :=Initialization.NewtonStep(
    //        x1[nR + 1:nF],
    //        b - C*invD*a,
    //        I - C*invD*B,
    //        nX);

    //      assert(sum(x2[nR+1:nF]) > 0, "x2[nR+1:nF] is zero at iteration " + String(index));
    //     x2 := P_to_orig_*Initialization.NewtonStep(P_to_id_*x1,f,J,nF);
    //         x2 := Initialization.NewtonStep(x1,f,J,nF);

        x1 :=x2;
         index :=index + 1;

        f :=Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps or index>1000 then
          solutionfound:=true;
        end if;
      end while;
      x2 :=P_to_orig_*x2;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      for i in 1:nF loop
      x_[i] :=max(1e-30,Xfull[i]);
      end for;
    //    x_ :=x2;//Xfull;

       success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_J.mat",
              "J",
              J);

       success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_f.mat",
              "f",
              {f});
      success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_J_.mat",
              "J_",
              J_);

    //  success :=Modelica.Utilities.Streams.writeRealMatrix(
    //         "U:/OutputX_pTXred_x1.mat",
    //         "x1",
    //         {x1});
      //
      // success :=Modelica.Utilities.Streams.writeRealMatrix(
      //        "U:/OutputX_pTXred_x1_.mat",
      //        "x1_",
      //        {x1});

    //   success :=Modelica.Utilities.Streams.writeRealMatrix(
    //     "U:/OutputX_pTXred_B_.mat",
    //     "B_",
    //     B);
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_D_.mat",
         "D_",
         D);
    //
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_n.mat",
         "n",
         {cat(1,pivots,{n_p,n_n})});
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_P.mat",
         "P",
         P);

    end X_pTXred_v4;

    function X_pTXred_v5
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-8;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B_;
      Real[nR] b_;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      MassFraction[nX] Xred_;

      Integer index;
      Boolean success;
      Real[nF,nX] lambda_mass1;

      Real[nR,nR] D;
      Real[nR] Dvector;
      Real[nR] invDvector;
      Real[nR,nX] B;
      Real[nX,nR] C;
      Real[nX,nX] I;
      Real[nR] a;
      Real[nX] b;
      Real[nR,nR] invD;
      Boolean[nR] I_p;
      Integer n_p;
      Integer n_n;
      Integer[nR] pivots;
      Integer[nR,nR] P;
      Real[nF,nF] J_;
      Real[nF] f_;
    algorithm
        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B_[r,firstnonzero[r]] :=1;
           b_[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xred,B_,b_);
         assert(sum(Xfull_in)> 0,"Test");
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;
    //        x1[1:nF] :=cat(1,{0.001*i for i in 1:nF-1},{1-0.001*sum(i for i in 1:nF-1)});//fill(1/nF, nF);
        x1[1:nF] :=fill(1/nF, nF);

         lambda_mass1 :=P_to_orig_*lambda_mass;
         // Xred_ :=Modelica.Math.Matrices.solve(transpose(lambda_mass),Xred);//
         Xred_:=transpose(lambda_mass_)*Xfull_in;

      // Newton algorithm
      solutionfound :=false;
      x1 :=P_to_id_*x1;
       index :=1;
      while not solutionfound loop
    //     f :=Initialization.calc_f.SLE_Tp_v1(x1,Xred,T,p);
    //     J :=Initialization.calc_J.SLE_Tp_v1(x1,T,p);
    //     x2 := Initialization.NewtonStep(x1,f,J,nF);
    //     assert(sum(x1)>0, "x1 is zero at iteration " + String(index));
        f :=-Initialization.calc_f.SLE_Tp(x1,Xred_,T,p);
        J :=Initialization.calc_J.SLE_Tp(x1,T,p);

        // check convergence
        if Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) < eps or index>100 then
          solutionfound:=true;
        end if;
         a :=f[1:nR];
         b :=f[nR + 1:nF];
         D :=J[1:nR, 1:nR];
         B :=J[1:nR, nR + 1:nF];
         C :=J[nR + 1:nF, 1:nR];
         I :=identity(nX);
         I_p :={false for i in 1:nR};//{abs(D[i, i]) >= Modelica.Math.Vectors.norm(C[:, i],Modelica.Constants.inf) for i in 1:nR};//

         n_p :=Modelica.Math.BooleanVectors.countTrue(I_p);
         pivots[1:n_p] :=Modelica.Math.BooleanVectors.index(I_p);
         pivots[n_p+1:nR] :=Modelica.Math.BooleanVectors.index(not I_p);
         n_n :=nR - n_p;
        P :=zeros(nR, nR);
         for i in 1:nR loop
           P[pivots[i],i] :=1;
         end for;
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_a.mat",
    //           "a",
    //           {a});
          a :=transpose(P)*a;
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_a_.mat",
    //           "a_",
    //           {a});
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_D.mat",
    //           "D",
    //           D);

    //      D :=P*D;
         Dvector :=D*fill(1, nR);
         Dvector :=Dvector[pivots];
         D :=diagonal(Dvector);
         C :=C*P;
    //       success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_B.mat",
    //           "B",
    //           B);
         B :=transpose(P)*B;
    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_x1.mat",
    //           "x1",
    //           {x1});
         x1[1:nR] :=transpose(P)*x1[1:nR];
         x2[1:nR] :=transpose(P)*x2[1:nR];
    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //           "U:/OutputX_pTXred_x1_.mat",
    //           "x1_",
    //           {x1});

    //      Dvector :=D*fill(1, nR);
          invDvector :=ones(nR) ./ Dvector;
         invD :=diagonal(invDvector);//fill(invDvector,nR);//diagonal(ones(nR)./Dvector);//{{if i == j then 1/D[i, j] else 0 for i in 1:nR} for j in 1:nR};
    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //          "U:/OutputX_pTXred_invD.mat",
    //          "invD",
    //          invD);

         J_[1:nR-n_p,1:nR-n_p] :=D[n_p + 1:nR, n_p + 1:nR];
         J_[1:nR-n_p,nR-n_p+1:n_n+nX] :=B[n_p + 1:nR, :];
         J_[nR-n_p+1:n_n+nX,1:nR-n_p] :=C[:, n_p + 1:nR];
            J_[nR-n_p+1:n_n+nX,nR-n_p+1:n_n+nX] :=I - C[:, 1:n_p]*invD[1:n_p, 1:n_p]*B[1:n_p, :]; // CORRECT

               for i in 1:nX loop
                 for j in 1:nX loop
                   J_[nR-n_p+i,nR-n_p+j] :=I[i, j] - C[i, 1:n_p]*(invDvector[1:n_p].*B[1:n_p,j]);
                 end for;
               end for;

    //       J_[nR-n_p+1:n_n+nX,nR-n_p+1:n_n+nX] :=identity(nX) - C[:, 1:n_p]*invDvector[1:n_p]*B[1:n_p, :];

         f_[1:nR-n_p] :=a[n_p + 1:nR];
    //   f_[nR-n_p+1:nF-n_p] :=b - C[:, 1:n_p]*invD[1:n_p, 1:n_p]*a[1:n_p]; // CORRECT
               for i in 1:nX loop
                  f_[nR-n_p+i] :=b[i] - C[i, 1:n_p]*(invDvector[1:n_p].*a[1:n_p]);//invD[1:n_p,1:n_p]*a[1:n_p];
    //              f_[nR-n_p+i] :=b[i] - C[i, 1:n_p]*invD[1:n_p,1:n_p]*a[1:n_p];
               end for;
    //       f_[nR-n_p+1:nF-n_p] :=b - C[:, 1:n_p]*(a[1:n_p].*invDvector[1:n_p]);//*(a[1:n_p].*invDvector[1:n_p]);

        // creating Permutation matrix and do pivoting, then solve subproblem with Newton

    //        x2[1:n_p] :=invD[1:n_p, 1:n_p]*(a[1:n_p] - B[1:n_p]*x1[nR + 1:nF]);
    //        while Modelica.Math.Vectors.norm(x2[n_p+1:nF]-x1[n_p+1:nF],Modelica.Constants.inf) > eps loop
           x2[n_p+1:nF] :=Initialization.NewtonStep(
            x1[n_p + 1:nF],
            f_[1:nF - n_p],
            J_[1:nF - n_p, 1:nF - n_p],
            nF - n_p);
    //          x1[n_p+1:nF] :=x2[n_p + 1:nF];
    //        end while;
    //         for i in 1:n_p loop
    //           x2[i] :=invD[i, i]*(a[i] - B[i, :]*x2[nR + 1:nF]);
    //         end for;

              x2[1:n_p] :=invD[1:n_p, 1:n_p]*(a[1:n_p] - B[1:n_p,:]*x2[nR + 1:nF]);//CORRECT

    //         x2[1:n_p] :=invDvector[1:n_p].*(a[1:n_p] - B[1:n_p]*x2[nR + 1:nF]);

           x2[1:nR] :=P*x2[1:nR];//*transpose(P);
    //        x1[1:nR] :=P*x1[1:nR];//*transpose(P);

    //      success :=Modelica.Utilities.Streams.writeRealMatrix(
    //          "U:/OutputX_pTXred_x2.mat",
    //          "x2",
    //          {x2});

    //      x2[1:nR] :=invD*(a - B*x1[nR + 1:nF]);
    //        x2[nR+1:nF] :=Initialization.NewtonStep(
    //        x1[nR + 1:nF],
    //        b - C*invD*a,
    //        I - C*invD*B,
    //        nX);

    //      assert(sum(x2[nR+1:nF]) > 0, "x2[nR+1:nF] is zero at iteration " + String(index));
    //     x2 := P_to_orig_*Initialization.NewtonStep(P_to_id_*x1,f,J,nF);
    //         x2 := Initialization.NewtonStep(x1,f,J,nF);

        x1 :=x2;
         index :=index + 1;
      end while;
      x2 :=P_to_orig_*x2;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      for i in 1:nF loop
      x_[i] :=max(1e-20,Xfull[i]);
      end for;
       x_ :=Xfull;

       success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_J.mat",
              "J",
              J);

       success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_f.mat",
              "f",
              {f});

       success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_f_.mat",
              "f_",
              {f_});
      success :=Modelica.Utilities.Streams.writeRealMatrix(
              "U:/OutputX_pTXred_J_.mat",
              "J_",
              J_);

    //  success :=Modelica.Utilities.Streams.writeRealMatrix(
    //         "U:/OutputX_pTXred_x1.mat",
    //         "x1",
    //         {x1});
      //
      // success :=Modelica.Utilities.Streams.writeRealMatrix(
      //        "U:/OutputX_pTXred_x1_.mat",
      //        "x1_",
      //        {x1});

       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_invD.mat",
         "invD",
         invD);
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_D_.mat",
         "D_",
         D);
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_Xred.mat",
         "Xred",
         {Xred_});
    //    success :=Modelica.Utilities.Streams.writeRealMatrix(
    //      "U:/OutputX_pTXred_b.mat",
    //      "b",
    //      {b});
    //
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_n.mat",
         "n",
         {cat(1,pivots,{n_p,n_n},b)});

       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_P_nu.mat",
         "P",
         P_to_id_);

       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_lambda_id.mat",
         "lambda_id",
         lambda_id);

       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_lambda_.mat",
         "lambda_",
         lambda_);
       success :=Modelica.Utilities.Streams.writeRealMatrix(
         "U:/OutputX_pTXred_P.mat",
         "P",
         P);

    end X_pTXred_v5;

    function X_pTXred_v6
      "SLE initialization according with reduced mass fraction as closed-system constraint. Solution algorithm with Newton solver according to Leal et al. (2016)"
      input SI.Pressure p;
      input SI.Temperature T;
      input MassFraction[nX] Xred;
      input Real[nF] Xinit = zeros(nF);
      output Real[nF] x_;

    protected
      parameter Real init = 1e-10;
      parameter Real eps = 1e-8;

      MassFraction[nX] Xrednonzero;
      MassFraction[nX] Xred_;

      Real[nF,nF] J;
      Real[nF] f;
      Real[nF] x1;
      Real[nF] x2;

      Boolean solutionfound = false;

      parameter Integer[nR] firstnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({nu[r,i] <> 0 for i in 1:nF}) for r in 1:nR};
      parameter Integer[nR] secondnonzero = {Modelica.Math.BooleanVectors.firstTrueIndex({if i>firstnonzero[r] then nu[r,i] <> 0 else false for i in 1:nF}) for r in 1:nR};

      Real[nR,nF] B;
      Real[nR] b;

      SI.MassFraction[nF] Xfull_in;
      SI.AmountOfSubstance[nF] nfull_in;
      SI.Mass[nF] mfull_in;

      SI.AmountOfSubstance[nF] ninit;
      SI.Mass[nF] minit;
      Real scale;

      SI.MoleFraction[nF] Yfull;
      SI.MassFraction[nF] Xfull;

      Integer index;
    algorithm
      //ensure no species amount to be zero for numerical issues
      for i in 1:nX loop
        Xrednonzero[i] :=max(Modelica.Constants.eps, Xred[i]);
      end for;

        // find initial guess x1 for Newton solver
         for r in 1:nR loop
           B[r,firstnonzero[r]] :=1;
           b[r] :=0;
         end for;
         Xfull_in :=Media.Common.equalityLeastSquares(transpose(lambda_mass),Xrednonzero,B,b);
         // nfull_in with positive elements only to fulfill transpose(lambda_mass)*mfull_in = Xred;
         if sum(Xfull_in) > 0 then
           mfull_in :=Xfull_in;
           nfull_in :=mfull_in ./ MMX;
           for i in 1:nF loop
             x1[i] :=max(init, nfull_in[i]);
           end for;
         elseif sum(Xinit)> 0 then
         // ninit from initial guess provided as input (Xinit)
           scale :=sum(transpose(lambda_mass_)*Xinit);
           minit :=Xinit*scale;
           ninit :=minit./MMX;
           for i in 1:nF loop
             x1[i] :=max(init, ninit[i]);
           end for;
         else
      // else generic initial guess
        x1[1:nF] :=fill(1/nF, nF);
         end if;

         Xred_:=transpose(lambda_mass_)*Xfull_in;

      // Newton algorithm
      solutionfound :=false;
      x1 :=P_to_id_*x1;
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
      x2 :=P_to_orig_*x2;
      Yfull :=x2[1:nF]/sum(x2[1:nF]);
      Xfull :=moleToMassFractions(Yfull, MMX);

      x_ :=Xfull;
    //   for i in 1:nF loop
    //   x_[i] :=max(1e-20,Xfull[i]);
    //   end for;

    end X_pTXred_v6;
  end MixtureLiquidAlgorithm;
end Common;
