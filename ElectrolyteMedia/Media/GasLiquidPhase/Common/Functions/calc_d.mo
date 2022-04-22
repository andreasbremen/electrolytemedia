within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions;
function calc_d
  "calculates density of gas, solute and solvent mixture at T and p"
  input SI.Temperature T;
  input SI.Density dg;
  input SI.MassFraction[nGfunfun+nLfunfun] X;

  output SI.Density d;
protected
  SI.SpecificVolume v;
algorithm

  v :=calc_v(
    T,
    dg,
    X);
  d :=1/v;

  annotation(smoothOrder=5);
end calc_d;
