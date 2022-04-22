within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_x_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] x_ijtautau;
protected
  Real A_phitautau=DebyeHueckel.Reduced.Additional.d2g_dtau2.calc_d2Aln_dtau2(T,p);
  Real I = calc_I(X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      x_ijtautau [i,j]:= 6*abs(datafun[i].z*datafun[j].z)*A_phitautau*I^0.5;
    end for;
  end for;
end calc_x_ijtautau;
