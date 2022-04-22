within ElectrolyteMedia.Media.SolidPhase.Common.Functions;
function calc_MMX "Vector of molar masses from each species"
  output SI.MolarMass[nSfun] MMX;

algorithm
  MMX[:] :=datafun[:].MM;
end calc_MMX;
