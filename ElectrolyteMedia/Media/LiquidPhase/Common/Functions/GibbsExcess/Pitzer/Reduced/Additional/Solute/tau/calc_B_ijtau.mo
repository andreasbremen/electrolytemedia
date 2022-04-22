within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_B_ijtau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] B_ijtau;
protected
  Real[nLifun,nLifun] beta0tau=calc_beta0_ijtau(T);
  Real[nLifun,nLifun] beta1tau=calc_beta1_ijtau(T);
  Real[nLifun,nLifun] beta2tau=calc_beta2_ijtau(T);
  Real[nLifun,nLifun] g_alpha1=GibbsExcess.Pitzer.Solute.calc_g(
                                             1, X);
  Real[nLifun,nLifun] g_alpha2=GibbsExcess.Pitzer.Solute.calc_g(
                                             2, X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if abs(datafun[i].z) > 0 and abs(datafun[j].z) > 0 then
        B_ijtau [i,j] := beta0tau[i,j] + beta1tau[i,j] * g_alpha1[i,j] + beta2tau[i,j] * g_alpha2[i,j];
      end if;
    end for;
  end for;
end calc_B_ijtau;
