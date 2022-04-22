within ElectrolyteMedia.Media.GasLiquidPhase.Common.MixtureLiquid;
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

  Real[1] f = {1};
  Real[1,1] J;
  SI.Temperature T_;
  SI.Temperature[1] Ttemp;

  SI.MassFraction[nG] Xg;
  SI.Density dg;

algorithm
  //calculate start values for inner-outer algorithm
  T :=exp((log(Tmin) + log(Tmax))/2);
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
    f := Initialization.calc_f.Entropy(T,dg,s,X);
    J := Initialization.calc_J.Entropy(T,dg,X);
    Ttemp:=Initialization.NewtonStep(
      {T},
      f,
      J,
      1);
    assert(sum(Ttemp)>0,"Error in Newton step for enthalpy calculation",AssertionLevel.error);
    T_ :=T;
    T :=Ttemp[1];

    if T <Tmin then
      T :=Tmin;
      T_ :=exp((log(Tmin) + log(Tmax))/2);
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

end XT_psXred;
