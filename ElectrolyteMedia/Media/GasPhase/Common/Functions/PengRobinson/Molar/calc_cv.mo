within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_cv "Calculates molar heat capacity of gas phase at constant volume with Peng Robinson"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output SI.MolarHeatCapacity cv;
protected
  Real p_dT(unit="Pa/K") = calc_der_p_T(
    T,
    d,
    y_i);
  Real v_dT(unit="m3/(K.mol)") = Molar.calc_v_dT_p(T,d,y_i);
  SI.MolarHeatCapacity cp = Molar.calc_cp(T,d,y_i);
algorithm
  cv :=cp - T*v_dT*p_dT;

end calc_cv;
