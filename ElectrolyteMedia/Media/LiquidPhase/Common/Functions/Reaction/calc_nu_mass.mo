within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.Reaction;
function calc_nu_mass "Mass based stoichiometry"
  input Real[:,:] nu;

  output Real[size(nu,1),size(nu,2)] nu_mass;

algorithm
  for i in 1:size(nu,2) loop
    if i < size(nu,2) then
      nu_mass[:,i] :=nu[:, i]*datafun[i].MM;
    else
      nu_mass[:,i] :=nu[:, i]*IF97.MH2O;
    end if;
  end for;

end calc_nu_mass;
