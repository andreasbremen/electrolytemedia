within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_beta "Calculates isobaric expansion coefficient of gas phase with corresponding EOS "
  input SI.Temperature T;
  input SI.Density d;
  input SI.MassFraction[:] X;

  output Real beta(unit="1/K");

algorithm
  if GasModelfun == Media.Common.Types.GasModel.Ideal then
    beta :=IdealGas.calc_beta(T);
  elseif GasModelfun == Media.Common.Types.GasModel.PengRobinson then
    beta :=PengRobinson.calc_beta(T,d,X);
  end if;

end calc_beta;
