within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_kappa "Calculates isothermal compressibility of gas phase with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output Real kappa(unit="1/Pa");
protected
  SI.Pressure p;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    p :=IdealGas.calc_p(T, d, X);
    kappa :=IdealGas.calc_kappa(p);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    kappa :=PengRobinson.calc_kappa(T,d,X);
  end if;

end calc_kappa;
