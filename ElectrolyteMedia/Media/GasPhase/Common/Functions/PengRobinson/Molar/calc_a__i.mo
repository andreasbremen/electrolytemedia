within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a__i
  "Calculates product of van der Waals attraction and acentric factor in Peng Robinson EOS"

  input SI.Temperature T;
  output Real[nGfun] a__i(unit="N.m4/mol2");
protected
  Real[nGfun] alpha_i(unit="") = Molar.calc_alpha_i(T);
  Real[nGfun] a_i(unit="N.m4/mol2") = Molar.calc_a_i();
algorithm
  for i in 1:nGfun loop
    a__i[i] :=a_i[i]*alpha_i[i];
  end for;

end calc_a__i;
