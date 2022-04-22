within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pipi;
function calc_x_ijpipi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] x_ijpipi;
protected
  Real A_phipipi=GibbsExcess.DebyeHueckel.Reduced.Additional.d2g_dpi2.calc_d2Aln_dpi2(T,p);
  Real I = calc_I(X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      x_ijpipi [i,j]:= 6*abs(datafun[i].z*datafun[j].z)*A_phipipi*I^0.5;
    end for;
  end for;
end calc_x_ijpipi;
