within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solute;
function calc_DH_contribution
  "Calculates first term in Bromley equation A_gamma*I^0.5/(1+I^0.5)"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real DebyeHuckel;
protected
  Real A_gamma=DebyeHueckel.calc_A_log( T,p);
  Real I=calc_I(X);
algorithm
  DebyeHuckel :=A_gamma*I^0.5/(1 + I^0.5);
end calc_DH_contribution;
