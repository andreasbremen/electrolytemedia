within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_B_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] B_ijtautau;
protected
  Real[nLifun,nLifun] beta0tautau=calc_beta0_ijtautau(
                                                T);
  Real[nLifun,nLifun] beta1tautau=calc_beta1_ijtautau(
                                                T);
  Real[nLifun,nLifun] beta2tautau=calc_beta2_ijtautau(
                                                T);
  Real[nLifun,nLifun] g_alpha1=GibbsExcess.Pitzer.Solute.calc_g(
                                             1, X);
  Real[nLifun,nLifun] g_alpha2=GibbsExcess.Pitzer.Solute.calc_g(
                                             2, X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if abs(datafun[i].z) > 0 and abs(datafun[j].z) > 0 then
        B_ijtautau [i,j] := beta0tautau[i,j] + beta1tautau[i,j] * g_alpha1[i,j] + beta2tautau[i,j] * g_alpha2[i,j];
      end if;
    end for;
  end for;
end calc_B_ijtautau;
