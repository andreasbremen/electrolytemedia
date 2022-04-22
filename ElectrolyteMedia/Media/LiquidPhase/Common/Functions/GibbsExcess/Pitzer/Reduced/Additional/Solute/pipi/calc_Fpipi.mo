within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pipi;
function calc_Fpipi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real Fpipi;

protected
  Real[nLifun] m_i=calc_mfromX(X);
  Real[nLifun,nLifun] der_E_theta_ijpipi=calc_der_E_theta_ijpipi(
      T,
      p,
      X);
  Real f_gammapipi=calc_f_gammapipi(
      T,
      p,
      X);

algorithm
  Fpipi := f_gammapipi;

  for i in 1:nLifun loop
    for j in i:nLifun loop
      Fpipi :=Fpipi + m_i[i]*m_i[j]*(der_E_theta_ijpipi[i, j]);
    end for;
  end for;

end calc_Fpipi;
