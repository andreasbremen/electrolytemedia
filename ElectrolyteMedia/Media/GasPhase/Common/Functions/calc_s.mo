within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_s "Calculates specific entropy of gas phase with corresponding EOS"
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output SI.SpecificEntropy s;

protected
  SI.Pressure p;

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    p :=IdealGas.calc_p(T, d, X);
    s :=IdealGas.calc_s(
      T,
      p,
      X);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    s :=PengRobinson.calc_s(T,d,X);
  end if;

end calc_s;
