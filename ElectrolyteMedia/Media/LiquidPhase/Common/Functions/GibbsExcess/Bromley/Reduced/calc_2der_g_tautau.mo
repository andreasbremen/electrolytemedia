within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced;
function calc_2der_g_tautau
  "Calculates 2nd derivative w.r.t. tau of excess reduced Gibbs free energy at infinite dilution of aqueous species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun] gtautau;

protected
  Real[nLfun-1] mol = calc_mfromX(X);
  Real[nLfun - 1] loggamma=Solute.calc_loggamma(T,p,X);
  Real loggamma_s = Solvent.calc_loggamma(T,p,X);
  Real[nLfun - 1] loggammatau=Additional.Solute.tau.calc_loggammatau(T,p,X);
  Real loggammatau_s = Additional.Solvent.tau.calc_loggammatau(T,p,X);
  Real[nLfun - 1] loggammatautau= Additional.Solute.tautau.calc_loggammatautau(T,p,X);
  Real loggammatautau_s = Additional.Solvent.tautau.calc_loggammatautau(T,p,X);

algorithm

  for i in 1:nLfun - 1 loop
    if mol[i] > 0 then
      gtautau[i] := -(10^loggamma[i]*mol[i])^(-2)*(mol[i]*log(10)*10^loggamma[i]*loggammatau[i])^2 + (10^loggamma[i]*mol[i])^(-1)
      *mol[i]*log(10)*(log(10)*10^loggamma[i]*loggammatau[i]^2 + 10^loggamma[i]*loggammatautau[i]);
    end if;
  end for;

  gtautau[nLfun] := -(10^loggamma_s*1/IF97.MH2O)^(-2)*(1/IF97.MH2O*log(10)*10^loggamma_s*loggammatau_s)^2 + (10^loggamma_s*1/IF97.MH2O)^(-1)
    *1/IF97.MH2O*log(10)*(log(10)*10^loggamma_s*loggammatau_s^2 + 10^loggamma_s*loggammatautau_s);

end calc_2der_g_tautau;
