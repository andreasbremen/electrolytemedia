within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_phi
  "Calculates fugacity coefficients of gas species with Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output Real[nGfun] phi;
protected
  Real[nGfun] lnphi=calc_lnphi(T,d,X);
algorithm
  for i in 1:nGfun loop
    phi[i] := exp(lnphi[i]);
  end for;

annotation(smoothOrder=5);
end calc_phi;
