within ElectrolyteMedia.Media.LiquidPhase.Common.MixtureLiquid;
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

  MoleFraction[nX] Yred;

  Real[1] f_T;
  Real[1,1] J_T;
  SI.Temperature T1;
  SI.Temperature T2;
  SI.Temperature[1] Tt_temp;

  Real[nF] f_X;
  Real[nF,nF] J_X;
  MoleFraction[nF] Y1;
  MoleFraction[nF] Y2;
algorithm
  //calculate start values for inner-outer algorithm
  T1 :=exp((log(Tmin) + log(Tmax))/2);
  X :=X_pTXred(p,T1,Xred);
  Y1 :=Functions.calc_Y(X);
  Yred :=transpose(lambda)*Y1;

  //inner (equilibrium) - outer (enthalpy) algorithm
  while abs(T2/T1-1) > eps loop
    f_T := Initialization.calc_f.Enthalpy(T1,p,h,X);
    J_T := Initialization.calc_J.Enthalpy(T1,p,X);
    Tt_temp:= Initialization.NewtonStep({T1},f_T,J_T,1);
    assert(sum(Tt_temp)>0,"Error in Newton step for enthalpy calculation",AssertionLevel.error);
    T2 :=T1;
    T1 :=Tt_temp[1];

    X :=X_pTXred(p,T1,Xred,X);
  end while;

  T :=T1;
end XT_phXred;
