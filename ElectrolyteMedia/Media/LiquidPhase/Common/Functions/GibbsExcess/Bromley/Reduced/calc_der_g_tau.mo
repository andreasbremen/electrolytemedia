within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced;
function calc_der_g_tau
  "Calculates derivative w.r.t. tau of excess reduced Gibbs free energy at infinite dilution of aqueous species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun] gtau;

protected
  Real[nLfun] ai=calc_ai(T,p,X);
  Real[nLfun-1] mol = calc_mfromX(X);
  Real[nLfun - 1] loggamma=Solute.calc_loggamma(T,p,X);
  Real loggamma_s = Solvent.calc_loggamma(T,p,X);
  Real[nLfun - 1] loggammatau=Additional.Solute.tau.calc_loggammatau(T,p,X);
  Real loggammatau_s = Additional.Solvent.tau.calc_loggammatau(T,p,X);

algorithm

  for i in 1:nLfun-1 loop
    if ai[i] > 0 then
      gtau[i] := 1/ai[i] * mol[i]*log(10)*10^loggamma[i]*loggammatau[i];
    end if;
  end for;

  gtau[nLfun] := 1/ai[nLfun]*1/IF97.MH2O*log(10)*10^loggamma_s*loggammatau_s;

end calc_der_g_tau;
