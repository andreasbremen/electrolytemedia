within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a_ij_d2T
  "Calculates second T derivative of product of van der Waals attraction and acentric factor in Peng Robinson EOS"

  input SI.Temperature T;

  output Real[nGfun,nGfun] a_ij_dT(unit="N.m4/(mol2.K2)");
protected
  Real[nGfun,nGfun] k_ij(unit="") = Molar.calc_k_ij(T);
  Real[nGfun,nGfun] k_ij_dT(unit="1/K")=interactionfun.k2_ij;
  Real[nGfun] a__i(unit="N.m4/mol2") = Molar.calc_a__i(T);
  Real[nGfun] a__i_dT(unit="N.m4/(mol2.K)") = Molar.calc_a__i_dT(T);
  Real[nGfun] a__i_d2T(unit="N.m4/(mol2.K2)") = Molar.calc_a__i_d2T(T);

  Real b;
  Real c;
  Real b_dT;
  Real c_dT;
  Real b_d2T;
  Real c_d2T;
algorithm
  for i in 1:nGfun loop
    for j in 1:nGfun loop
      b :=a__i[i];
      c :=a__i[j];
      b_dT :=a__i_dT[i];
      c_dT :=a__i_dT[j];
      b_d2T :=a__i_d2T[i];
      c_d2T :=a__i_d2T[j];
      a_ij_dT[i,j] := -2*k_ij_dT[i,j]*((c)^0.5*b_dT/(2*(b)^0.5) + (b)^0.5*c_dT/(2*(c)^0.5)) + (1-k_ij[i,j])*(b_dT*c_dT/(2*(b*c)^0.5) + (c)^0.5*(b_d2T/(2*(b)^0.5) - (b_dT)^2/(4*(b)^(3/2))) + (b)^0.5*(c_d2T/(2*(c)^0.5) - (c_dT)^2/(4*(c)^(3/2))));
    end for;
  end for;

end calc_a_ij_d2T;
