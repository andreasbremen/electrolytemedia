within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a_d2T
  "Calculates second T derivative of product of van der Waals attraction and acentric factor in Peng Robinson EOS"
  input SI.Temperature T;
  input SI.MoleFraction[nGfun] y_i;

  output Real a_dT(unit="N.m4/(mol2.K2)");
protected
  Real[nGfun,nGfun] a_ij_d2T(unit="N.m4/(mol2.K2)") = Molar.calc_a_ij_d2T(T);
algorithm
  for i in 1:nGfun loop
    for j in 1:nGfun loop
      a_dT := a_dT + y_i[i]*y_i[j]*a_ij_d2T[i, j];
    end for;
  end for;
end calc_a_d2T;
