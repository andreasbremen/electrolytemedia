within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_kappa_i
  "Species specific constant correlated against acentric factor"
  output Real[nGfun] kappa_i(unit="");
algorithm
  for i in 1:nGfun loop
    kappa_i[i] :=0.37464 + 1.54226 * datafun[i].w - 0.26992 * (datafun[i].w)^2;
  end for;
end calc_kappa_i;
