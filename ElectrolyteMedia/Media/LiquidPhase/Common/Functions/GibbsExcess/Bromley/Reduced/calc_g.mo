within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced;
function calc_g "Calculates excess reduced Gibbs free energy at infinite dilution of aqueous species"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun] g_red;

protected
  Real[nLfun] g=Bromley.calc_g(T,p,X);

algorithm

  for i in 1:nLfun loop
  g_red[i] := g[i]/(Modelica.Constants.R*T);
  end for;

end calc_g;
