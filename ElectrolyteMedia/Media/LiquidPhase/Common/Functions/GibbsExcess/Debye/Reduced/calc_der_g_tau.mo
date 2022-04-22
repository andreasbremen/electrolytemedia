within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Debye.Reduced;
function calc_der_g_tau
  "Calculates derivative w.r.t. tau of excess reduced Gibbs free energy at infinite dilution of aqueous species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction X[nLfun];
  output Real gtau[nLfun];

protected
  Real[nLfun] a_i = calc_activity(T,p,X);
  Real[nLfun-1] mol_i = calc_mfromX(X);
  Real[nLfun] loggamma = calc_loggamma(T,p,X);
  Real Atau = DebyeHueckel.Reduced.Additional.dg_dtau.calc_dAlog_dtau(T,p);
  Real I = calc_I(X);
  Real mol_s = sum(mol_i);

algorithm

   for i in 1:nLfun-1 loop
     if a_i[i] > 0 then
       gtau[i] := 1/(a_i[i]) * mol_i[i]*log(10)*10^loggamma[i]*(-Atau)*datafun[i].z^2*sqrt(I);
     end if;
   end for;
   gtau[nLfun] := 1/(a_i[nLfun]) * 1/MixtureLiquid.MH2O*log(10)*10^loggamma[nLfun]*log10(exp(1))*IF97.MH2O*mol_s*1/3*sqrt(I)*Atau;

end calc_der_g_tau;
