within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_Y_j

  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1] Y_j;
protected
  Real mol_i[nLfun-1] =  calc_mfromX(X);
  Real I_i[nLfun-1];
  Real I_a=calc_I_a(X);
algorithm
  for i in 1:nLfun-1 loop
    if datafun[i].z < 0 then
      I_i[i] :=0.5*mol_i[i]*datafun[i].z^2;
    else
      I_i[i] :=0;
    end if;
  end for;

  Y_j :=I_i./(I_a+1e-25);

end calc_Y_j;
