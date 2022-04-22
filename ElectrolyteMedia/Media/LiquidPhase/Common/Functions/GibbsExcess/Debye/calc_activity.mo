within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Debye;
function calc_activity
  "Calculates activity in molality base of aqueous species from Debye limiting law"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Real[nLfun] a;
protected
  Real[nLfun] gamma;
  Real[nLfun] m;
algorithm

  m[1:nLfun-1] :=calc_mfromX(X);
  m[nLfun] := 1/MixtureLiquid.MH2O;
  gamma :=calc_gamma(
          T,
          p,
          X);
  a:=gamma .* m;
  annotation(smoothOrder=20);
end calc_activity;
