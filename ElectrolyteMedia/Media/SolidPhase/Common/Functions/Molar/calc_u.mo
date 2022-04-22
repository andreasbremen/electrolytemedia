within ElectrolyteMedia.Media.SolidPhase.Common.Functions.Molar;
function calc_u "Calculates molar internal energy of mineral phases"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MoleFraction[nSfun] x;

  output SI.MolarEnergy u;
protected
  SI.MolarEnergy u_i[nSfun]=Molar.calc_u_i(T,p);
algorithm
  u :=u_i*x;
  annotation(smoothOrder=2);
end calc_u;
