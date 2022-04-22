within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced.Sym;
function GibbsDerivs "Returns record of values of dimensionless Gibbs function and its derivatives at given temperature and pressure"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Modelica.Media.Common.GibbsDerivs gderivs;
protected
  Real tau=calc_tau(T);
  Real pi=calc_pi(p);
  Real g=calc_g(T,X);
  Real gpi=calc_der_g_pi();
  Real gpipi=calc_2der_g_pipi();
  Real gtau=calc_der_g_tau();
  Real gtautau=calc_2der_g_tautau();
  Real gtaupi=calc_2der_g_taupi();
  SI.SpecificEntropy R = ElectrolyteMedia.Media.LiquidPhase.Common.Functions.calc_R( X);

algorithm

  gderivs.T := T;
  gderivs.p := p;
  gderivs.tau := tau;
  gderivs.pi := pi;
  gderivs.g := g;
  gderivs.gpi := gpi;
  gderivs.gpipi := gpipi;
  gderivs.gtau := gtau;
  gderivs.gtautau := gtautau;
  gderivs.gtaupi := gtaupi;
  gderivs.R := R;

end GibbsDerivs;
