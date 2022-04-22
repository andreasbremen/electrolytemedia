within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pi;
function calc_x_ijpi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real [nLifun, nLifun] x_ijpi;
protected
  Real A_phipi=GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dpi.calc_dAln_dpi(T,p);
  Real I = calc_I(X);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      x_ijpi [i,j]:= 6*abs(datafun[i].z*datafun[j].z)*A_phipi*I^0.5;
    end for;
  end for;
end calc_x_ijpi;
