within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_alpha_i
  "Calculates van der Waals attraction of each species in Peng Robinson EOS"
  input SI.Temperature T;
  output Real[nGfun] alpha_i(unit="");
protected
  Real[nGfun] kappa_i=Molar.calc_kappa_i();
algorithm
  for i in 1:nGfun loop
    alpha_i[i] :=(1 + kappa_i[i]*(1 - sqrt(T/datafun[i].T_c)))^2;
  end for;
end calc_alpha_i;
