within ElectrolyteMedia.Media.Common.Reaction;
function calc_lambda_id
  "calculates lambda with nX*nX identity matrix"
  input Real[:,:] nu_id;
  output Real[size(nu_id,2),size(nu_id,2)-size(nu_id,1)] lambda_id;
algorithm
  lambda_id :=transpose(cat(2,-transpose(nu_id[:, size(nu_id,1) + 1:size(nu_id,2)]),identity(size(nu_id,2)-size(nu_id,1))));
end calc_lambda_id;
