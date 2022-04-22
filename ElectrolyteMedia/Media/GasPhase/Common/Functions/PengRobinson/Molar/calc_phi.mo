within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_phi
  "Calculates fugacity coefficients of gas species with Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y;

  output Real[nGfun] phi;
protected
  Real[nGfun] lnphi=Molar.calc_lnphi(T,d,y);
algorithm
  for i in 1:nGfun loop
    phi[i] := exp(lnphi[i]);
  end for;
end calc_phi;
