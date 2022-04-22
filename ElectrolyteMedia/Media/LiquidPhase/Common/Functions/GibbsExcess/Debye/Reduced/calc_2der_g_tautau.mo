within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Debye.Reduced;
function calc_2der_g_tautau
  "Calculates 2nd derivative w.r.t. tau of excess reduced Gibbs free energy at infinite dilution of aqueous species"

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
  Real Atautau = DebyeHueckel.Reduced.Additional.d2g_dtau2.calc_d2Alog_dtau2(T,p);
  Real I = calc_I(X);

algorithm

   for i in 1:nLfun-1 loop
     if a_i[i] > 0 then
       gtautau[i] := 1/a_i[i]*mol_i[i]*log(10)*(-10^loggamma[i]*datafun[i].z^2*sqrt(I)*Atautau+log(10)*10^loggamma[i]*(-Atau*datafun[i].z^2*sqrt(I))^2) -1/(a_i[i])^2 * (mol_i[i]*log(10)*10^loggamma[i]*(-Atau)*datafun[i].z^2*sqrt(I))^2;
     end if;
   end for;
   gtautau[nLfun] := 1/a_i[nLfun]*1/MixtureLiquid.MH2O*log(10)*(10^loggamma[nLfun]*log10(exp(1))*IF97.MH2O*mol_s*1/3*sqrt(I)*Atautau+log(10)*10^loggamma[nLfun]*(log10(exp(1))*IF97.MH2O*mol_s*1/3*sqrt(I)*Atau)^2)-1/a_i[nLfun]^2*(1/MixtureLiquid.MH2O*log(10)*10^loggamma[nLfun]*log10(exp(1))*IF97.MH2O*mol_s*1/3*sqrt(I)*Atau)^2;

end calc_2der_g_tautau;
