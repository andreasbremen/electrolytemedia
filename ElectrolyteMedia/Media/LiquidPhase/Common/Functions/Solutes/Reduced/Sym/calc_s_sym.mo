within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Solutes.Reduced.Sym;
function calc_s_sym "Calculates reduced entropy for asymmetrical mixing"
  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output SI.SpecificEntropy s;
protected
  SI.SpecificEntropy R = ElectrolyteMedia.Media.LiquidPhase.Common.Functions.calc_R( X);
algorithm

      s :=R*(calc_tau(T)*calc_der_g_tau() - calc_g(T, X));

end calc_s_sym;
