within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solvent;
function calc_B_Phi_ij
  "Calculates cation i - anion j interaction B_Phi_ij"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] B_Phi_ij;
protected
  Real [nLifun,nLifun] beta0 = Solute.calc_beta0_ij(T);
  Real [nLifun,nLifun] beta1 = Solute.calc_beta1_ij(T);
  Real [nLifun,nLifun] beta2 = Solute.calc_beta2_ij(T);
  Real alpha1[nLifun, nLifun] = Solute.calc_alpha1();
  Real alpha2[nLifun, nLifun] = Solute.calc_alpha2();
  Real I = calc_I(X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if abs(datafun[i].z) > 0 and abs(datafun[j].z) > 0 then
        B_Phi_ij [i,j] := beta0[i,j] + beta1 [i,j] * exp(-alpha1[i,j]*I^0.5) + beta2[i,j] * exp(-alpha2[i,j]*I^0.5);
      end if;
    end for;
  end for;
end calc_B_Phi_ij;
