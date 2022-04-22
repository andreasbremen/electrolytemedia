within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_k_ij
  "Calculates temperature dependent interactionfun parameter in Peng Robinson EOS"

  input SI.Temperature T;

  output Real[nGfun,nGfun] k_ij(unit="");
algorithm

 k_ij := interactionfun.k1_ij + interactionfun.k2_ij*T;
end calc_k_ij;
