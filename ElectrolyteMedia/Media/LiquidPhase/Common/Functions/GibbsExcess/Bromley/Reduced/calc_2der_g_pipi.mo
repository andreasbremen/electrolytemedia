within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced;
function calc_2der_g_pipi
  "Calculates 2nd derivative w.r.t. pi of excess reduced Gibbs free energy at infinite dilution of aqueous species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun] gpipi;

protected
  Real[nLfun-1] mol = calc_mfromX(X);
  Real[nLfun - 1] loggamma=Solute.calc_loggamma(T,p,X);
  Real loggamma_s=Solvent.calc_loggamma(T,p,X);
  Real[nLfun - 1] loggammapi=Additional.Solute.pi.calc_loggammapi(T,p,X);
  Real loggammapi_s = Additional.Solvent.pi.calc_loggammapi(T,p,X);
  Real[nLfun - 1] loggammapipi=Additional.Solute.pipi.calc_loggammapipi(T,p,X);
  Real loggammapipi_s=Additional.Solvent.pipi.calc_loggammapipi(T,p,X);

algorithm

 for i in 1:nLfun-1 loop
   if mol[i] > 0 then
     gpipi[i] := -(10^loggamma[i]*mol[i])^(-2)*(mol[i]*log(10)*10^loggamma[i]*loggammapi[i])^2 + (10^loggamma[i]*mol[i])^(-1)
      *mol[i]*log(10)*(log(10)*10^loggamma[i]*loggammapi[i]^2 + 10^loggamma[i]*loggammapipi[i]);
   end if;
 end for;

 gpipi[nLfun] := -(10^loggamma_s*1/IF97.MH2O)^(-2)*(1/IF97.MH2O*log(10)*10^loggamma_s*loggammapi_s)^2 + (10^loggamma_s*1/IF97.MH2O)^(-1)
    *1/IF97.MH2O*log(10)*(log(10)*10^loggamma_s*loggammapi_s^2 + 10^loggamma_s*loggammapipi_s);

end calc_2der_g_pipi;
