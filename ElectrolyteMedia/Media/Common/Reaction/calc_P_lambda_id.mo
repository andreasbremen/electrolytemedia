within ElectrolyteMedia.Media.Common.Reaction;
function calc_P_lambda_id
  "Calculates a stoichiometry matrix with pivoting so that lambdaout =(I|lambda_var) with permutation matrix Pout for reordering species elements of lambda_ordered = lambdaout*Pout"
  input Real[:,:] lambda;
  output Integer[size(lambda,1),size(lambda,1)] Pout;
protected
  Real[size(lambda,1),size(lambda,2)] lambdaout;
algorithm
  (lambdaout,Pout) :=calc_lambda_id(lambda);

end calc_P_lambda_id;
