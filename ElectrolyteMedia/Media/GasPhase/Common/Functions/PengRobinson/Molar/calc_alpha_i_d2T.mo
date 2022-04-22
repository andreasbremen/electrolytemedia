within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_alpha_i_d2T
  "Calculates second T derivative of van der Waals attraction of each species in Peng Robinson EOS"
  input SI.Temperature T;
  output Real[nGfun] alpha_i_d2T(unit="1/K2");
protected
  Real[nGfun] kappa_i=Molar.calc_kappa_i();
algorithm
  alpha_i_d2T :=kappa_i.*(kappa_i+ones(nGfun)).*sqrt(T./datafun.T_c)/(2*T^2);
end calc_alpha_i_d2T;
