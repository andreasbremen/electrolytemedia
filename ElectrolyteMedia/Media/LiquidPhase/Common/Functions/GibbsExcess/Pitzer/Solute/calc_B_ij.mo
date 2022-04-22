within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_B_ij "Calculates cation i - anion j interaction B_ij"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] B_ij;
protected
  Real [nLifun,nLifun] beta0 = calc_beta0_ij(T);
  Real [nLifun,nLifun] beta1 = calc_beta1_ij(T);
  Real [nLifun,nLifun] beta2 = calc_beta2_ij(T);
  Real[nLifun,nLifun] g_alpha1=Solute.calc_g(1, X);
  Real[nLifun,nLifun] g_alpha2=Solute.calc_g(2, X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if abs(datafun[i].z) > 0 and abs(datafun[j].z) > 0 then
        B_ij [i,j] := beta0[i,j] + beta1 [i,j] * g_alpha1[i,j] + beta2[i,j] * g_alpha2[i,j];
      end if;
    end for;
  end for;
end calc_B_ij;
