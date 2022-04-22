within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pitau;
function calc_Fpitau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real Fpitau;

protected
  Real[nLifun] m_i=calc_mfromX(X);
  Real[nLifun,nLifun] der_E_theta_ijpitau=calc_der_E_theta_ijpitau(
      T,
      p,
      X);
  Real f_gammapitau=calc_f_gammapitau(
      T,
      p,
      X);

algorithm
  Fpitau := f_gammapitau;

  for i in 1:nLifun loop
    for j in i:nLifun loop
      Fpitau :=Fpitau + m_i[i]*m_i[j]*(der_E_theta_ijpitau[i, j]);
    end for;
  end for;

end calc_Fpitau;
