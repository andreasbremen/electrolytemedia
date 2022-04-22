within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Debye.Reduced;
function calc_2der_g_taupi
  "Calculates 2nd derivative w.r.t. tau,pi of excess reduced Gibbs free energy at infinite dilution of aqueous species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real gtautau[nLfun];

protected
  Real[nLfun] a_i = calc_activity(T,p,X);
  Real[nLfun-1] mol_i = calc_mfromX(X);
  Real mol_s = sum(mol_i);
  Real[nLfun] loggamma = calc_loggamma(T,p,X);
  Real Atau = DebyeHueckel.Reduced.Additional.dg_dtau.calc_dAlog_dtau(T,p);
  Real Api = DebyeHueckel.Reduced.Additional.dg_dpi.calc_dAlog_dpi(T,p);
  Real Ataupi = DebyeHueckel.Reduced.Additional.d2g_dtaudpi.calc_d2Alog_dtaudpi(T,p);
  Real I = calc_I(X);

algorithm

   for i in 1:nLfun-1 loop
     if a_i[i] > 0 then
       gtautau[i] := 1/a_i[i]*mol_i[i]*log(10)*(-10^loggamma[i]*datafun[i].z^2*sqrt(I)*Ataupi+log(10)*10^loggamma[i]*(datafun[i].z^2*sqrt(I))^2*Atau*Api) -1/(a_i[i])^2 * (mol_i[i]*log(10)*10^loggamma[i]*datafun[i].z^2*sqrt(I))^2*Atau*Api;
     end if;
   end for;
   gtautau[nLfun] := 1/a_i[nLfun]*1/MixtureLiquid.MH2O*log(10)*(10^loggamma[nLfun]*log10(exp(1))*IF97.MH2O*mol_s*1/3*sqrt(I)*Ataupi+log(10)*10^loggamma[nLfun]*(log10(exp(1))*IF97.MH2O*mol_s*1/3*sqrt(I))^2*Atau*Api)-1/a_i[nLfun]^2*(1/MixtureLiquid.MH2O*log(10)*10^loggamma[nLfun]*log10(exp(1))*IF97.MH2O*mol_s*1/3*sqrt(I))^2*Atau*Api;

end calc_2der_g_taupi;
