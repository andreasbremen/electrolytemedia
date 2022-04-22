within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_der_v_T
  "Calculates derivative w.r.t. T at const. p of specific volume of gas mixture with Peng Robinson EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[nGfun] X;

  output Real v_dT(unit="m3/(kg.K)");
protected
  SI.MoleFraction[nGfun] y=calc_Y(X);
algorithm
  v_dT := Molar.calc_v_dT_p(T,d,y)/calc_MM(X);

end calc_der_v_T;
