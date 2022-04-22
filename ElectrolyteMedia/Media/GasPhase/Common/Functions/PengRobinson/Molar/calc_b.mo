within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_b "Calculates gas minimal volume b in Peng Robinson EOS"

  input SI.MoleFraction[nGfun] y_i;
  output SI.MolarVolume b;
protected
  SI.MolarVolume[nGfun] b_i=Molar.calc_b_i();
algorithm
  b := y_i * b_i;
end calc_b;

function calc_B
  "Calculates derivative of compressibility factor of minimal gas volume in Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output Real B(unit="");
protected
  SI.Pressure p=Molar.calc_p(T,d,y_i);
  SI.MolarVolume b=Molar.calc_b(y_i);
algorithm
  B := b*p/(Modelica.Constants.R*T);
end calc_B;
