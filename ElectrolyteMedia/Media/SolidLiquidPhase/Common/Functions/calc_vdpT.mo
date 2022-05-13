within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions;
function calc_vdpT
  "calculates pressure derivative of specific volume of solid mixture at T and p"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nSfunfun+nLfunfun] X;

  output Real vdpT;

protected
  SI.MassFraction[nSfunfun] Xs = calc_Xs(X);
  SI.MassFraction[nLfunfun] Xl = calc_Xl(X);
  SI.DerDensityByPressure dsdpT;
  SI.DerDensityByPressure dldpT;
  Real vsdpT;
  Real vldpT;
  SI.Density ds;
  SI.Density dl;
algorithm
  ds :=SolidFunctions.calc_d(T,p,Xs);
  dl :=LiquidFunctions.calc_d(T,p,Xl);

  dsdpT :=SolidPhase.Common.Functions.calc_der_d_p(
    T,
    p,
    Xs);
  dldpT :=LiquidPhase.Common.Functions.calc_der_d_p(
    T,
    p,
    Xl);

  vsdpT :=-1/ds^2*dsdpT;
  vldpT :=-1/dl^2*dldpT;

  vdpT :=sum(X[1:nSfunfun])*vsdpT + sum(X[1 + nSfunfun:nSfunfun + nLfunfun])*vldpT;

end calc_vdpT;
