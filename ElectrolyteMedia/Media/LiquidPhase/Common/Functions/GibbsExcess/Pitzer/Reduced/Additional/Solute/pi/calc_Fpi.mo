within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pi;
function calc_Fpi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real Fpi;

protected
  Real[nLifun] m_i=calc_mfromX(X);
  Real[nLifun,nLifun] der_E_theta_ijpi=calc_der_E_theta_ijpi(
      T,
      p,
      X);
  Real f_gammapi=calc_f_gammapi(
      T,
      p,
      X);

algorithm
  Fpi := f_gammapi;

  for i in 1:nLifun loop
    for j in i:nLifun loop
      Fpi :=Fpi + m_i[i]*m_i[j]*(der_E_theta_ijpi[i, j]);
    end for;
  end for;

end calc_Fpi;
