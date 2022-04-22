within ElectrolyteMedia.Media.Common.Reaction;
function calc_P_nu_id
  "Calculates a stoichiometry matrix with pivoting so that nuout =(I|nu_var) with permutation matrix Pout for reordering species elements of nu_ordered = nuout*Pout"
  input Real[:,:] nu;
  output Integer[size(nu,2),size(nu,2)] Pout;
protected
  Real[size(nu,1),size(nu,2)] nuout;
algorithm
  (nuout,Pout) :=calc_nu_id(nu);

end calc_P_nu_id;
