within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_ddpT
  "calculates pressure derivative of density of solid mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output SI.DerDensityByPressure ddpT;

protected
  SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
  SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
  SI.DerDensityByPressure dsdpT;
  SI.DerDensityByPressure dldpT;
  Real vsdpT;
  Real vldpT;
  Real vdpT;
  SI.Density ds;
  SI.Density dl;
  SI.Density d;
algorithm
  ds :=SolidFunctions.calc_d(T,p,Xs);
  dl :=LiquidFunctions.calc_d(T,p,Xl);

  dsdpT :=SolidFunctions.calc_der_d_p(
    T,
    p,
    Xs);
  dldpT :=LiquidFunctions.calc_der_d_p(
    T,
    p,
    Xl);

  vsdpT :=-1/ds^2*dsdpT;
  vldpT :=-1/dl^2*dldpT;

  vdpT :=sum(X[1:nSfunfun])*vsdpT + sum(X[1 + nSfunfun:nSfunfun + nLfunfun])*vldpT;

  d :=calc_d(T, p, X);

  ddpT :=-vdpT*d^2;

end calc_ddpT;
