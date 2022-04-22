within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_der_v_p
  "calculates pressure derivative of specific volume of solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Density dg;
  input SI.MassFraction[nGfunfun+nLfunfun] X;

  output Real vdpT;

protected
  SI.MassFraction[nGfunfun] Xg = X[1:nGfunfun]/sum(X[1:nGfunfun]);
  SI.Pressure p = GasFunctions.calc_p(T,dg,Xg);
  SI.MassFraction[nLfunfun] Xl = X[1+nGfunfun:nGfunfun+nLfunfun]/sum(X[1+nGfunfun:nGfunfun+nLfunfun]);
  SI.DerDensityByPressure dgdpT;
  SI.DerDensityByPressure dldpT;
  Real vgdpT;
  Real vldpT;
//   SI.Density dg;
  SI.Density dl;
algorithm
//   dg :=GasFunctions.calc_d(T,p,Xg);
  dl :=LiquidFunctions.calc_d(T,p,Xl);

  dgdpT :=GasFunctions.calc_der_d_p(
    T,
    dg,
    Xg);
  dldpT :=LiquidFunctions.calc_der_d_p(
    T,
    p,
    Xl);

  vgdpT :=-1/dg^2*dgdpT;
  vldpT :=-1/dl^2*dldpT;

  vdpT :=sum(X[1:nGfunfun])*vgdpT + sum(X[1 + nGfunfun:nGfunfun + nLfunfun])*vldpT;

end calc_der_v_p;
