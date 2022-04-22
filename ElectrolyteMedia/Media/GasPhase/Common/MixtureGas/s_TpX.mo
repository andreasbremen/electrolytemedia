within ElectrolyteMedia.Media.GasPhase.Common.MixtureGas;
function s_TpX "Return specific entropy"
  extends Modelica.Icons.Function;
  input Temperature T "Temperature";
  input AbsolutePressure p;
  input MassFraction X[:]=reference_X
    "Independent Mass fractions of gas mixture";
  output SI.SpecificEntropy s "Specific entropy at temperature T";

protected
  Density d = d_TpX(T,p,X);

algorithm
  s :=Functions.calc_s(
      T,
      d,
      X);
  annotation (Inline=false, smoothOrder=2);
end s_TpX;
