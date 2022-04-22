within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_F "Pitzer, p. 89"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real F;
protected
  Real [nLifun] m_i = calc_mfromX(X);
  Real [nLifun,nLifun] beta0 = calc_beta0_ij(T);
  Real [nLifun,nLifun] beta1 = calc_beta1_ij(T);
  Real [nLifun,nLifun] beta2 = calc_beta2_ij(T);
  Real[nLifun,nLifun] der_E_theta_ij=Solute.calc_der_E_theta_ij(
      T,
      p,
      X);
  Real[nLifun,nLifun] der_B_ij=Solute.calc_der_B_ij(T, X);
  Real f_gamma=Solute.calc_f_gamma(
      T,
      p,
      X);
algorithm
  F := f_gamma;

  for i in 1:nLifun loop
    for j in i:nLifun loop
      F :=F + m_i[i]*m_i[j]*(der_B_ij[i, j] + der_E_theta_ij[i, j]);
    end for;
  end for;

end calc_F;
