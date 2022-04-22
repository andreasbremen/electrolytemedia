within ElectrolyteMedia.Media.GasPhase.Common.MixtureGas;
function g_TpX "Return specific enthalpy"
  extends Modelica.Icons.Function;
  input SI.Temperature T "Temperature";
  input SI.Pressure p "Pressure";
  input MassFraction X[:]=reference_X "Independent Mass fractions of gas mixture";
  output SI.SpecificGibbsFreeEnergy g "Specific Gibbs free energy at temperature T";
protected
  Density d = d_TpX(T,p,X);
algorithm
  g :=Functions.calc_g(
      T,
      d,
      X);
  annotation(Inline=false,smoothOrder=2);
end g_TpX;
