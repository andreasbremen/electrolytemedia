within ElectrolyteMedia.Media.SolidLiquidPhase.Common.Functions.Reaction;
function calc_nu_mass
  input Real[:,:] nu;

  output Real[size(nu,1),size(nu,2)] nu_mass;
protected
  SI.MolarMass[size(nu,2)] MMX = calc_MMX();
algorithm
  for i in 1:size(nu,1) loop
    nu_mass[i,:] :=nu[i, :] .* MMX[:];
  end for;

//   if GasModelfunfun == Media.Common.Types.GasModel.Ideal then
//     for i in 1:nSfunfun loop
//       nu_mass[:,i] :=nu[:, i] * dataIGfunfun[i].MM;
//     end for;
//   elseif GasModelfunfun == Media.Common.Types.GasModel.PengRobinson then
//     for i in 1:nSfunfun loop
//       nu_mass[:,i] :=nu[:,i] * dataPRfunfun[i].MM;
//     end for;
//   end if;
//   for i in 1:nSfunfun:nSfunfun+nLfunfun-1 loop
//     nu_mass[:,i] :=nu[:,i] * dataLfunfun[i].MM;
//   end for;
//   nu_mass[:,nSfunfun+nLfunfun] :=nu[:, nSfunfun + nLfunfun] * IF97.MH2O;
end calc_nu_mass;
