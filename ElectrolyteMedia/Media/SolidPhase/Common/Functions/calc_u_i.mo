within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_u_i
  "Calculates specific internal energy of mineral species"

  input SI.Temperature T;
  input SI.Pressure p;

  output SI.SpecificEnthalpy u_i[nSfun];

protected
  SI.MolarEnthalpy[nSfun] u_i_m = Molar.calc_u_i(T,p);

algorithm
    u_i :=u_i_m./datafun[:].MM;
  annotation(smoothOrder=5);
end calc_u_i;
