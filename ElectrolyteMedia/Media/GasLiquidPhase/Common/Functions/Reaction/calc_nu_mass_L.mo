within ElectrolyteMedia.Media.GasLiquidPhase.Common.Functions.Reaction;
function calc_nu_mass_L
  input Real[:,:] nu_L;

  output Real[size(nu_L,1),size(nu_L,2)] nu_mass_L;
protected
  SI.MolarMass[:] MMX_L = LiquidFunctions.calc_MMX();
algorithm
  if size(nu_L,2)>0 then
  for i in 1:size(nu_L,1) loop
    nu_mass_L[i,:] :=nu_L[i, :] .* MMX_L[:];//MMX[1+nGfunfun:nGfunfun:nLfunfun];
  end for;
  end if;

//   if GasModelfunfun == Media.Common.Types.GasModel.Ideal then
//     for i in 1:nGfunfun loop
//       nu_mass[:,i] :=nu[:, i] * dataIGfunfun[i].MM;
//     end for;
//   elseif GasModelfunfun == Media.Common.Types.GasModel.PengRobinson then
//     for i in 1:nGfunfun loop
//       nu_mass[:,i] :=nu[:,i] * dataPRfunfun[i].MM;
//     end for;
//   end if;
//   for i in 1:nGfunfun:nGfunfun+nLfunfun-1 loop
//     nu_mass[:,i] :=nu[:,i] * dataLfunfun[i].MM;
//   end for;
//   nu_mass[:,nGfunfun+nLfunfun] :=nu[:, nGfunfun + nLfunfun] * IF97.MH2O;
end calc_nu_mass_L;
