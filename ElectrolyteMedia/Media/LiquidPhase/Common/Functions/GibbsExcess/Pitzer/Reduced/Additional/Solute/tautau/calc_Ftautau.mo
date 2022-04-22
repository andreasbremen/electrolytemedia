within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_Ftautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real Ftautau;
protected
  Real [nLifun] m_i = calc_mfromX(X);
  Real[nLifun,nLifun] der_E_theta_ijtautau=calc_der_E_theta_ijtautau(
      T,
      p,
      X);
  Real[nLifun,nLifun] der_B_ijtautau=calc_der_B_ijtautau(
                                                   T, X);
  Real f_gammatautau=calc_f_gammatautau(
          T,
          p,
          X);
algorithm
  Ftautau := f_gammatautau;

  for i in 1:nLifun loop
    for j in i:nLifun loop
      Ftautau :=Ftautau + m_i[i]*m_i[j]*(der_B_ijtautau[i, j] + der_E_theta_ijtautau[i, j]);
    end for;
  end for;

end calc_Ftautau;
