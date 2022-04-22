within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.tau;
function calc_B_Phi_ijtau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] B_Phi_ijtau;
protected
  Real [nLifun,nLifun] beta0tau = GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_beta0_ijtau(T);
  Real [nLifun,nLifun] beta1tau = GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_beta1_ijtau(T);
  Real [nLifun,nLifun] beta2tau = GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_beta2_ijtau(T);
  Real alpha1[nLifun, nLifun] = GibbsExcess.Pitzer.Solute.calc_alpha1();
  Real alpha2[nLifun, nLifun] = GibbsExcess.Pitzer.Solute.calc_alpha2();
  Real I = calc_I(X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if abs(datafun[i].z) > 0 and abs(datafun[j].z) > 0 then
        B_Phi_ijtau [i,j] := beta0tau[i,j] + beta1tau[i,j] * exp(-alpha1[i,j]*I^0.5) + beta2tau[i,j] * exp(-alpha2[i,j]*I^0.5);
      end if;
    end for;
  end for;
end calc_B_Phi_ijtau;
