within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_alpha_i_dT
  "Calculates T derivative of van der Waals attraction of each species in Peng Robinson EOS"
  input SI.Temperature T;
  output Real[nGfun] alpha_i_dT(unit="1/K");
protected
  Real[nGfun] kappa_i=Molar.calc_kappa_i();
algorithm
   for i in 1:nGfun loop
     alpha_i_dT[i] :=-(kappa_i[i]*sqrt(T/datafun[i].T_c)*(1-kappa_i[i]*(sqrt(T/datafun[i].T_c)-1)))/T;
   end for;
end calc_alpha_i_dT;
