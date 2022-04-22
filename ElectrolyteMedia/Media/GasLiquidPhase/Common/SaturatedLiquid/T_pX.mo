within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid;
function T_pX
  "Calculates boiling temperature of liquid solution at pressure p and overall composition X"
  input SI.Pressure p;
  input SI.MassFraction[nG+nL] X;
  output SI.Temperature T;
protected
  SI.MassFraction[nG] Xg = X[1:nG]/sum(X[1:nG]);
  SI.MassFraction[nL] Xl = X[1+nG:nG+nL]/sum(X[1+nG:nG+nL]);
  SI.MoleFraction[nL] Yl = Functions.LiquidFunctions.calc_Y(Xl);
  SI.Temperature T1;
  SI.Temperature T2;
  SI.SpecificEnergy gi[nG+nL];
  SI.MolarEnergy gim[nG+nL];
  SI.MolarEnergy grw;
  Real Kw;
  Real Kw1;
  Real[nL] gamma_i;
  SI.Density dg;
  Real[nG] phi_i;
  Boolean solutionfound = false;
  Real eps = 1e-3;

algorithm
  T1 :=373.15;

  while not solutionfound loop
    gi :=Functions.calc_gi(T1, p);
    for i in 1:nG+nL loop
      gim[i] :=gi[i]*MMX[i];
    end for;
    grw :=gim[nG + nL] - gim[nG];
    Kw :=exp(-grw/Modelica.Constants.R/T1);
    gamma_i :=Functions.LiquidFunctions.calc_gamma(T1, p,Xl);
    dg :=dg_TpX(T1,p,Xg);
    phi_i :=Functions.GasFunctions.calc_phi(T1,dg,Xg);
    Kw1 :=gamma_i[nL]/MH2O*Yl[nL]*prefg/(phi_i[nG]*p);

    T2 :=T_Kwp(Kw1, p);

    if abs(T2-T1)/T1 < eps then
      solutionfound :=true;
    end if;

    T1 :=T2;
  end while;

  T :=T2;
end T_pX;
