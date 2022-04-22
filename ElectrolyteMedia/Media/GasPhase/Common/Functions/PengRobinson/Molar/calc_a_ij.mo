within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a_ij
  "Calculates product of van der Waals attraction and acentric factor in Peng Robinson EOS"

  input SI.Temperature T;

  output Real[nGfun,nGfun] a_ij(unit="N.m4/mol2");
protected
  Real[nGfun,nGfun] k_ij(unit="") = Molar.calc_k_ij(T);
  Real[nGfun] a__i(unit="N.m4/mol2") = Molar.calc_a__i(T);
algorithm
  for i in 1:nGfun loop
    for j in 1:nGfun loop
      a_ij[i,j] := (1-k_ij[i,j])*a__i[i]^0.5*a__i[j]^0.5;
    end for;
  end for;
end calc_a_ij;
