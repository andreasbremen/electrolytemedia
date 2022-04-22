within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.tautau;
function calc_B_Phi_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] B_Phi_ijtautau;
protected
  Real [nLifun,nLifun] beta0tautau = GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau.calc_beta0_ijtautau(T);
  Real [nLifun,nLifun] beta1tautau = GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau.calc_beta1_ijtautau(T);
  Real [nLifun,nLifun] beta2tautau = GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau.calc_beta2_ijtautau(T);
  Real alpha1[nLifun, nLifun] = GibbsExcess.Pitzer.Solute.calc_alpha1();
  Real alpha2[nLifun, nLifun] = GibbsExcess.Pitzer.Solute.calc_alpha2();
  Real I = calc_I(X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if abs(datafun[i].z) > 0 and abs(datafun[j].z) > 0 then
        B_Phi_ijtautau[i,j] := beta0tautau[i,j] + beta1tautau[i,j] * exp(-alpha1[i,j]*I^0.5) + beta2tautau[i,j] * exp(-alpha2[i,j]*I^0.5);
      end if;
    end for;
  end for;
end calc_B_Phi_ijtautau;
