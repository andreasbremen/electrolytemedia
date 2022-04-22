within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_log_a_ij "calculates decadic logarithm of water activity"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] log_a_i;
protected
  Real[nLfun - 1,nLfun - 1] ln_a_i=calc_ln_a_ij(T,p, X);
algorithm

  for i in 1:nLfun - 1 loop
    for j in 1:nLfun - 1 loop
      log_a_i[i, j] := 1/log(10)*ln_a_i[i, j];
    end for;
  end for;

end calc_log_a_ij;
