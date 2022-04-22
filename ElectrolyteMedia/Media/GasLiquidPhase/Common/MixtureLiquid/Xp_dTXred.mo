within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid;
function Xp_dTXred
  "Gas-liquid and liquid dissociation equilibrium calculated via Newton method with homotopy continuation method of step-by step adding solutes to initially pure aqueous solution"
  input Modelica.SIunits.Density d;
  input SI.Temperature T;
  input MassFraction[nX] Xred;
  output MassFraction[nG+nL] X;
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

  Real pscale = 1e-5;
  Boolean bool = false;

  MassFraction[nG] Xg;
  SI.Density dg;
algorithm
  //calculate start values for inner-outer algorithm
  p :=exp((log(pmin) + log(pmax))/2);
  X :=X_pTXred(
    p,
    T,
    Xred);

  Xg :=X[1:nG]/sum(X[1:nG]);
  dg :=dg_TpX(
    T,
    p,
    Xg);

  //inner (equilibrium) - outer (enthalpy) algorithm
  while Modelica.Math.Vectors.norm(f,Modelica.Constants.inf) > eps loop

    f := Initialization.calc_f.Density(T,dg,X,d);
    J := Initialization.calc_J.Density(T,dg,X);
    ptemp:=Initialization.NewtonStep(
      {p},
      f,
      J,
      1);
    assert(sum(ptemp)>0,"Error in Newton step for enthalpy calculation",AssertionLevel.error);
    p_ :=p;
    p :=ptemp[1];

    if p <pmin then
      p :=pmin;
      p_ :=exp((log(pmin) + log(pmax))/2);
    end if;

    //inner algorithm
    X :=X_pTXred(
      p,
      T,
      Xred,
      X);
    Xg :=X[1:nG]/sum(X[1:nG]);
    dg :=dg_TpX(
      T,
      p,
      Xg);
  end while;
end Xp_dTXred;
