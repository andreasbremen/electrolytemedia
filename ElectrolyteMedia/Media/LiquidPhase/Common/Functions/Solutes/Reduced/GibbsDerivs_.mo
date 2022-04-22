within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced;
function GibbsDerivs_
  "Gibbs derivatives with Href as reference for reduced Gibbs free energy calculation"
  input SI.Temperature T;
  input SI.Pressure p;
  output Modelica.Media.Common.GibbsDerivs[nLifun] gderivs;
protected
  Real tau = calc_tau(T);
  Real pi = calc_pi(p);
  Real[nLifun] g = calc_g_(T,p);
  Real[nLifun] gpi=calc_der_g_pi(T, p);
  Real[nLifun] gpipi=calc_2der_g_pipi(T, p);
  Real[nLifun] gtau=calc_der_g_tau(T, p);
  Real[nLifun] gtautau=calc_2der_g_tautau(T, p);
  Real[nLifun] gtaupi=calc_2der_g_taupi(T, p);
algorithm
  for i in 1:nLifun loop
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
    gderivs[i].R :=datafun[i].R;
  end for;

end GibbsDerivs_;
