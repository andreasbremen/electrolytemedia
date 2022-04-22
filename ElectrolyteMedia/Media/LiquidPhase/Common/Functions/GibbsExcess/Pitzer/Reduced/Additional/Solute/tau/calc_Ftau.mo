within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_Ftau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real Ftau;
protected
  Real [nLifun] m_i = calc_mfromX(X);
  Real[nLifun,nLifun] der_E_theta_ijtau=calc_der_E_theta_ijtau(
      T,
      p,
      X);
  Real[nLifun,nLifun] der_B_ijtau=calc_der_B_ijtau(T, X);
  Real f_gammatau=calc_f_gammatau(
          T,
          p,
          X);
algorithm
  Ftau := f_gammatau;

  for i in 1:nLifun loop
    for j in i:nLifun loop
      Ftau :=Ftau + m_i[i]*m_i[j]*(der_B_ijtau[i, j] + der_E_theta_ijtau[i, j]);
    end for;
  end for;

end calc_Ftau;
