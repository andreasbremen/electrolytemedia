within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pitau;
function calc_x_ijpitau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] x_ijpitau;
protected
  Real A_phipitau=GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dtaudpi.calc_d2Aln_dtaudpi(T,p);
  Real I = calc_I(X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      x_ijpitau [i,j]:= 6*abs(datafun[i].z*datafun[j].z)*A_phipitau*I^0.5;
    end for;
  end for;
end calc_x_ijpitau;
