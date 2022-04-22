within ElectrolyteMedia.Media.LiquidPhase.Common.Functions;
function calc_der_d_p
  "Calculates pressure derivative of density of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[:] X;
  output SI.DerDensityByPressure ddpT;

protected
  SI.Density[nLfun] di;
  SI.DerDensityByPressure[nLfun] ddpTi;
  Real[nLfun] vi;
  Real[nLfun] dvpTi;
algorithm
  di[1:nLfun-1] :=Solutes.calc_d_i(T, p);
  di[nLfun] :=IF97_R1_Tp.calc_rho(T, p);
  ddpTi[1:nLfun-1] :=Solutes.calc_der_d_p(T, p);
  ddpTi[nLfun] :=IF97_R1_Tp.calc_der_d_p(T, p);

  for i in 1:nLfun-1 loop
    if datafun[i].MM > 0.001008 then
      vi[i] := 1 / di[i];
      dvpTi[i] :=-ddpTi[i]/di[i]^2;
    else
      vi[i] :=0;
      dvpTi[i] :=0;
    end if;
  end for;
  vi[nLfun] :=1/di[nLfun];
  dvpTi[nLfun] :=-ddpTi[nLfun]/di[nLfun]^2;
  ddpT :=-sum(X[i]*dvpTi[i] for i in 1:nLfun)/(vi*X)^2;

end calc_der_d_p;
