within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid;
function X_phXred
  "Liquid dissociation equilibrium calculated via Newton method with homotopy continuation method of step-by step adding solutes to initially pure aqueous solution"
  input Modelica.SIunits.Pressure p;
  input SI.SpecificEnthalpy h;
  input MassFraction[nX] Xred;
  output Real[nF] X;

protected
  parameter Real eps = 1e-6;

  Real[1] f;
  Real[1,1] J;
  Real vfrac1;
  Real vfrac2;
  Real[1] vfractemp;

  SI.Temperature T;
  SI.MassFraction[nG] Xg;
  SI.Density dg;

algorithm
  //calculate start values for inner-outer algorithm
  vfrac1 :=0.1;
  X :=X_vfracpXred(vfrac1,p,Xred);
  T :=T_pX(p, X);

  Xg :=X[1:nG]/sum(X[1:nG]);
  dg :=dg_TpX(T,p,Xg);

  //inner (equilibrium) - outer (enthalpy) algorithm
  while Modelica.Math.Vectors.norm({vfrac2/vfrac1-1}) > eps loop
     f := Initialization.calc_f.Enthalpy(T,dg,h,X);
     J := Initialization.calc_J.Enthalpy_vfrac(T,dg,X);
     vfractemp:=Initialization.NewtonStep({vfrac1},f,J,1);
     assert(sum(vfractemp)>0,"Error in Newton step for vapor fraction calculation",AssertionLevel.error);
     vfrac2 :=vfrac1;
     vfrac1 :=vfractemp[1];

    //inner algorithm
    X :=X_vfracpXred(vfrac1,p,Xred,X);
    T :=T_pX(p, X);
    Xg :=X[1:nG]/sum(X[1:nG]);
    dg :=dg_TpX(T,p,Xg);
  end while;

end X_phXred;
