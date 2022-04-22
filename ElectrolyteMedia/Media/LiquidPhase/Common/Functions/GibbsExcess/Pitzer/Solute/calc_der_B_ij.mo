within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_der_B_ij
  "Calculates derivative of cation i - anion j interaction B_ij"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] B_ij;
protected
  Real [nLifun,nLifun] beta0 = calc_beta0_ij(T);
  Real [nLifun,nLifun] beta1 = calc_beta1_ij(T);
  Real [nLifun,nLifun] beta2 = calc_beta2_ij(T);
  Real I = calc_I(X);
  Real[nLifun,nLifun] der_g_alpha1=Solute.calc_der_g(1, X);
  Real[nLifun,nLifun] der_g_alpha2=Solute.calc_der_g(2, X);
algorithm
  if I > 0 then
    for i in 1:nLifun loop
      for j in 1:nLifun loop
        if abs(datafun[i].z) > 0 and abs(datafun[j].z) > 0 then
          B_ij [i,j] := (beta1 [i,j] * der_g_alpha1[i,j] + beta2[i,j] * der_g_alpha2[i,j]) / I;
        end if;
      end for;
    end for;
  end if;
end calc_der_B_ij;
