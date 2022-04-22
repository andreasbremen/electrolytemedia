within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a_ij_dT
  "Calculates T derivative of product of van der Waals attraction and acentric factor in Peng Robinson EOS"

  input SI.Temperature T;

  output Real[nGfun,nGfun]
    a_ij_dT(unit="N.m4/(mol2.K)");
protected
  Real[nGfun,nGfun] k_ij(unit="") = Molar.calc_k_ij(T);
  Real[nGfun,nGfun] k_ij_dT(unit="1/K")=interactionfun.k2_ij;
  Real[nGfun] a__i(unit="N.m4/mol2") = Molar.calc_a__i(T);
  Real[nGfun] a__i_dT(unit="N.m4/(mol2.K)")= Molar.calc_a__i_dT(T);
algorithm
  for i in 1:nGfun loop
    for j in 1:nGfun loop
      a_ij_dT[i,j] := (a__i[j]*(-(k_ij[i,j]-1)*a__i_dT[i]-2*a__i[i]*k_ij_dT[i,j])-a__i[i]*(k_ij[i,j]-1)*a__i_dT[j])/(2*a__i[i]^0.5*a__i[j]^0.5);
    end for;
  end for;

end calc_a_ij_dT;
