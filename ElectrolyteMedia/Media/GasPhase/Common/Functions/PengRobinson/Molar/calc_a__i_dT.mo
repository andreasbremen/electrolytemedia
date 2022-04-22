within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a__i_dT
  "Calculates T derivative of product of van der Waals attraction and acentric factor in Peng Robinson EOS"

  input SI.Temperature T;
  output Real[nGfun] a__i_dT(unit="N.m4/(mol2.K)");
protected
  Real[nGfun] alpha_i_dT(unit="1/K")= Molar.calc_alpha_i_dT(T);
  Real[nGfun] a_i(unit="N.m4/mol2") = Molar.calc_a_i();
algorithm
  for i in 1:nGfun loop
    a__i_dT[i] := a_i[i] * alpha_i_dT[i];
  end for;
end calc_a__i_dT;
