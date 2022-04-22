within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_x_ij
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] x_ij;
protected
  Real A_phi=DebyeHueckel.calc_A_ln(T,p);
  Real I = calc_I(X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      x_ij [i,j]:= 6*abs(datafun[i].z*datafun[j].z)*A_phi*I^0.5;
    end for;
  end for;
end calc_x_ij;
