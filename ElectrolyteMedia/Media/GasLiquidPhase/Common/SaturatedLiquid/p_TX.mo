within ElectrolyteMedia.Media.GasLiquidPhase.Common.SaturatedLiquid;
function p_TX
  "Calculates vapor pressure of liquid solution at temperature T and overall composition X"
  input SI.Temperature T;
  input SI.MassFraction[nG+nL] X;
  output SI.Pressure p;
protected
  SI.MassFraction[nG] Xg = X[1:nG]/sum(X[1:nG]);
  SI.MassFraction[nL] Xl = X[1+nG:nG+nL]/sum(X[1+nG:nG+nL]);
  SI.MoleFraction[nL] Yl = Functions.LiquidFunctions.calc_Y(Xl);
  SI.Pressure p1;
  SI.Pressure p2;
  SI.SpecificEnergy gi[nG+nL];
  SI.MolarEnergy gim[nG+nL];
  SI.MolarEnergy grw;
  Real Kw;
  Real[nL] gamma_i;
  SI.Density dg;
  Real[nG] phi_i;
  Boolean solutionfound = false;
  Real eps = 1e-9;

algorithm
  p1 :=2e5;
  while not solutionfound loop
    gi :=Functions.calc_gi(T, p1);
    for i in 1:nG+nL loop
      gim[i] :=gi[i]*MMX[i];
    end for;
    grw :=gim[nG + nL] - gim[nG];
    Kw :=exp(-grw/Modelica.Constants.R/T);
    gamma_i :=Functions.LiquidFunctions.calc_gamma(T, p1,Xl);
    dg :=dg_TpX(T,p1,Xg);
    phi_i :=Functions.GasFunctions.calc_phi(T,dg,Xg);
    p2 :=gamma_i[nL]/MH2O*Yl[nL]*prefg/(phi_i[nG]*Kw);
    if abs(p2-p1)/p1 < eps then
      solutionfound :=true;
    end if;
    p1 :=p2;
  end while;

  p :=p2;
end p_TX;
