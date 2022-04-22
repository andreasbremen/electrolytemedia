within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a__i_d2T
  "Calculatessecond T derivative of product of van der Waals attraction and acentric factor in Peng Robinson EOS"

  input SI.Temperature T;
  output Real[nGfun] a__i_dT(unit="N.m4/(mol2.K2)");
protected
  Real[nGfun] alpha_i_d2T(unit="1/K2")= Molar.calc_alpha_i_d2T(T);
  Real[nGfun] a_i(unit="N.m4/mol2") = Molar.calc_a_i();
algorithm
  for i in 1:nGfun loop
    a__i_dT[i] := a_i[i] * alpha_i_d2T[i];
  end for;
end calc_a__i_d2T;
