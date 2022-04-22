within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Reduced;
function GibbsDerivs "Returns record of values of dimensionless Gibbs function and its derivatives at given temperature and pressure"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Modelica.Media.Common.GibbsDerivs[nLfun] gderivs;
protected
  Real tau = calc_tau(T);
  Real pi = calc_pi(p);
  Real[nLfun] g = calc_g(T,p,X);
  Real[nLfun] gpi=calc_der_g_pi(
      T,
      p,
      X);
  Real[nLfun] gpipi=calc_2der_g_pipi(
      T,
      p,
      X);
  Real[nLfun] gtau=calc_der_g_tau(
      T,
      p,
      X);
  Real[nLfun] gtautau=calc_2der_g_tautau(
      T,
      p,
      X);
  Real[nLfun] gtaupi=calc_2der_g_taupi(
      T,
      p,
      X);
algorithm
  for i in 1:nLfun loop
    gderivs[i].T :=T;
    gderivs[i].p :=p;
    gderivs[i].tau :=tau;
    gderivs[i].pi :=pi;
    gderivs[i].g :=g[i];
    gderivs[i].gpi :=gpi[i];
    gderivs[i].gpipi :=gpipi[i];
    gderivs[i].gtau :=gtau[i];
    gderivs[i].gtautau :=gtautau[i];
    gderivs[i].gtaupi :=gtaupi[i];
    if i == nLfun then
      gderivs[i].R := MixtureLiquid.RH2O;
    else
      gderivs[i].R :=datafun[i].R;
    end if;
  end for;

end GibbsDerivs;
