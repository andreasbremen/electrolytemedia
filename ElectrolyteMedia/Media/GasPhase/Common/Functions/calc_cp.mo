within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_cp "Calculates specific heat capacity at constant pressure of gas phase with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output SI.SpecificHeatCapacityAtConstantPressure cp;

protected
  SI.Pressure p;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    p :=IdealGas.calc_p(T, d, X);
    cp :=IdealGas.calc_cp(T,p,X);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    cp :=PengRobinson.calc_cp(T,d,X);
  end if;

end calc_cp;
