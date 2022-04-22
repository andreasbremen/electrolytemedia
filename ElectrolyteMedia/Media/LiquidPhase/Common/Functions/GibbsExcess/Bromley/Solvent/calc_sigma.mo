within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_sigma

  input SI.MassFraction[nLfun] X "mole fraction/molality/concentration";

  output Real sigma;
protected
  Real I = calc_I(X);
  Real rho = 1; //given in Bromley, 1973
algorithm
  sigma :=3/(rho*I^1.5+1e-25)*(1 + rho*I^0.5 - 1/(1 + rho*I^0.5) - 2*log(
    1 + rho*I^0.5));
end calc_sigma;
