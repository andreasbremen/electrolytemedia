within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_f "calculates the fugacity of gas phase species"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output Real[nGfun] f;
protected
  SI.Pressure p = calc_p(T,d,X);
  SI.MoleFraction[nGfun] Y = calc_Y(X);
  Real[nGfun] phi = calc_phi(T,d,X);
algorithm

  f :=phi .* Y*p/1e5;
annotation(smoothOrder=5);
end calc_f;
