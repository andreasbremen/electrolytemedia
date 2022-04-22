within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_x_ijtau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] x_ijtau;
protected
  Real A_phitau=DebyeHueckel.Reduced.Additional.dg_dtau.calc_dAln_dtau(T,p);
  Real I = calc_I(X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      x_ijtau [i,j]:= 6*abs(datafun[i].z*datafun[j].z)*A_phitau*I^0.5;
    end for;
  end for;
end calc_x_ijtau;
