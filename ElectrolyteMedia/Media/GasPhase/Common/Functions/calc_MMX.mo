within ElectrolyteMedia.Media.GasPhase.Common.Functions;
function calc_MMX "Returns vector containing molar masses of all species"

  output SI.MolarMass[nGfun] MMX;

algorithm
    MMX :=datafun[:].MM;
annotation(smoothOrder=5);
end calc_MMX;
