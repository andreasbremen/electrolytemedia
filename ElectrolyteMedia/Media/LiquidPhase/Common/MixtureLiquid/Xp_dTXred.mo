within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid;
function Xp_dTXred
  "Liquid dissociation equilibrium calculated via Newton method with homotopy continuation method of step-by step adding solutes to initially pure aqueous solution"
  input SI.Density d;
  input SI.Temperature T;
  input MassFraction[nX] Xred;
  output Real[nF] X;
  output SI.Pressure p;

protected
  parameter Real eps = 1e-5;

  SI.Pressure pmin = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.psat(T);
  SI.Pressure pmax = 1000e5;

  MoleFraction[nX] Yred;

  Real[1] f_p;
  Real[1,1] J_p;
  SI.Pressure p1;
  SI.Pressure p2;
  SI.Temperature[1] p_temp;

  Real[nF] f_X;
  Real[nF,nF] J_X;
  MoleFraction[nF] Y1;
  MoleFraction[nF] Y2;
algorithm
  //calculate start values for inner-outer algorithm
  p1 :=exp((log(pmin) + log(pmax))/2);
  X :=X_pTXred(p1,T,Xred);
  Y1 :=Functions.calc_Y(X);
  Yred :=transpose(lambda)*Y1;

  //inner (equilibrium) - outer (pressure) algorithm
  while abs(p2/p1-1) > eps loop
    f_p := Initialization.calc_f.Density(T,p1, X,d);
    J_p := Initialization.calc_J.Density(T,p1,X);
    p_temp:= Initialization.NewtonStep({p1},f_p,J_p,1);
    assert(sum(p_temp)>0,"Error in Newton step for pressure calculation",AssertionLevel.error);
    p2 :=p1;
    p1 :=p_temp[1];

    X :=X_pTXred(p1,T,Xred,X);
  end while;

  p :=p1;
end Xp_dTXred;
