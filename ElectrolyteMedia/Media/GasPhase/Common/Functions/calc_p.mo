within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_p "Calculates pressure of gas phase with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output SI.Pressure p;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    p :=IdealGas.calc_p(T, d, X);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    p :=PengRobinson.calc_p(T,d,X);
  end if;

end calc_p;
