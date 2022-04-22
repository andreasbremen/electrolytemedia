within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid;
function X_dTXred
  "Gas-liquid and liquid dissociation equilibrium calculated via Newton method with homotopy continuation method of step-by step adding solutes to initially pure aqueous solution"
  input Modelica.SIunits.Density d;
  input SI.Temperature T;
  input MassFraction[nX] Xred;
  output MassFraction[nG+nL] X;

protected
  SI.Pressure p;
  parameter Real eps = 1e-5;

  SI.Pressure pmin = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.psat(T);
  SI.Pressure pmax = 1000e5;

  Real[1] f={1};
  Real[1,1] J;
  Real vfrac;
  SI.Pressure vfrac_;
  SI.Pressure[1] vfractemp;

  SI.SpecificVolume v = 1/d;

  Real pscale = 1e-5;
  Boolean bool = false;

  MassFraction[nG] Xg;
  SI.Density dg;
algorithm

  //calculate start values for inner-outer algorithm
  vfrac :=0.5;
  X :=X_vfracTXred(vfrac,T,Xred);
  p :=p_TX(T,X);

  Xg :=X[1:nG]/(sum(X[1:nG])+Modelica.Constants.eps);
  dg :=dg_TpX(T,p,Xg);

  //inner (equilibrium) - outer (enthalpy) algorithm
   while Modelica.Math.Vectors.norm({vfrac_/vfrac-1}) > eps loop
      f := Initialization.calc_f.Density(T,dg,X,d);
      J := Initialization.calc_J.Density_vfrac(T,dg,X);
      vfractemp:=Initialization.NewtonStep({vfrac},f,J,1);
      assert(sum(vfractemp)>0,"Error in Newton step for vapor fraction calculation",AssertionLevel.error);
      vfrac_ :=vfrac;
      vfrac :=vfractemp[1];

     //inner algorithm
     X :=X_vfracTXred(vfrac,T,Xred,X);
     p :=p_TX(T, X);
     Xg :=X[1:nG]/sum(X[1:nG]);
     dg :=dg_TpX(T,p,Xg);
 end while;
end X_dTXred;
