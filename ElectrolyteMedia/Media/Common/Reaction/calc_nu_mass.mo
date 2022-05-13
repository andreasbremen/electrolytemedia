within ElectrolyteMedia.Media.Common.Reaction;
function calc_nu_mass "Mass based stoichiometry"
  input Real[:,:] nu;
  input Real[size(nu,2)] MMX;

  output Real[size(nu,1),size(nu,2)] nu_mass;

algorithm
  for i in 1:size(nu,2) loop
    nu_mass[:, i] := nu[:, i]*MMX[i];
  end for;

end calc_nu_mass;
