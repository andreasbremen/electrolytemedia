within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a_i_sum
  "Calculates product of van der Waals attraction and acentric factor of species in Peng Robinson EOS"

  input SI.Temperature T;
  input SI.MoleFraction[nGfun] y_i;

  output Real[nGfun] a_i_sum(unit="N.m4/mol2");
protected
  Real[nGfun,nGfun] a_ij(unit="N.m4/mol2") = Molar.calc_a_ij(T);
algorithm
  a_i_sum := zeros(nGfun);
  for k in 1:nGfun loop
    for i in 1:nGfun loop
      a_i_sum[k] := a_i_sum[k] + y_i[i]*a_ij[i, k];
    end for;
  end for;
end calc_a_i_sum;
