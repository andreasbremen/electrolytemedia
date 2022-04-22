within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_cv "Calculates specific heat capacity at constant volume of gas phase with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output SI.SpecificHeatCapacityAtConstantVolume cv;

protected
  SI.Pressure p;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    p :=IdealGas.calc_p(T, d, X);
    cv :=IdealGas.calc_cv(T,p,X);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    cv :=PengRobinson.calc_cv(T,d,X);
  end if;

end calc_cv;
