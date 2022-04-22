within ElectrolyteMedia.Media.GasPhase.Common.MixtureGas;
function h_TpX "Return specific enthalpy"
  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input AbsolutePressure p;
  input MassFraction X[:]=reference_X
    "Independent Mass fractions of gas mixture";
  output SI.SpecificEnthalpy h "Specific enthalpy at temperature T";

protected
  Density d = d_TpX(T,p,X);

algorithm
  h :=Functions.calc_h(
      T,
      d,
      X);
  annotation (Inline=false, smoothOrder=2);
end h_TpX;
