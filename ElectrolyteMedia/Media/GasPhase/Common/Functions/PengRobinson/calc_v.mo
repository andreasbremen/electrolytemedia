within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson;
function calc_v
  "Calculates specific volume of gas phase with Peng Robinson"
  input SI.Density d;

  output SI.SpecificVolume v;

algorithm

  v:=1/d;

end calc_v;
