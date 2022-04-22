within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.IF97_R1_Tp.Reduced;
function GibbsDerivs "Returns record of values of dimensionless Gibbs function and its derivatives at given temperature and pressure"
  input SI.Temperature T;
  input SI.Pressure p;
  output Modelica.Media.Common.GibbsDerivs gderivs;
protected
  Real tau = calc_tau(T);
  Real pi = calc_pi(p);
  Real g = calc_g(T,p);
  Real gpi=calc_der_g_pi(T, p);
  Real gpipi=calc_2der_g_pipi(T, p);
  Real gtau=calc_der_g_tau(T, p);
  Real gtautau=calc_2der_g_tautau(T, p);
  Real gtaupi=calc_2der_g_taupi(T, p);
algorithm
  gderivs.T :=T;
  gderivs.p :=p;
  gderivs.tau :=tau;
  gderivs.pi :=pi;
  gderivs.g :=g;
  gderivs.gpi :=gpi;
  gderivs.gpipi :=gpipi;
  gderivs.gtau :=gtau;
  gderivs.gtautau :=gtautau;
  gderivs.gtaupi :=gtaupi;
  gderivs.R :=IF97.RH2O;

end GibbsDerivs;
