within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a_dT
  "Calculates T derivative of product of van der Waals attraction and acentric factor in Peng Robinson EOS"
  input SI.Temperature T;
  input SI.MoleFraction[nGfun] y_i;

  output Real a_dT(unit="N.m4/(mol2.K)");
protected
  Real[nGfun,nGfun] a_ij_dT(unit="N.m4/(mol2.K)")=Molar.calc_a_ij_dT(T);
  Real[nGfun,nGfun] a_ij_dT_(unit="N.m4/(mol2.K)");

algorithm
  for i in 1:nGfun loop
    for j in 1:nGfun loop
      a_ij_dT_[i,j] := y_i[i]*y_i[j]*a_ij_dT[i, j];
    end for;
  end for;
  a_dT :=sum(a_ij_dT_);
end calc_a_dT;
