within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_I_a

  input SI.MassFraction[nLfun] X;

  output Real I_a;

protected
  Real mol_i[nLfun-1] =  calc_mfromX(X);

algorithm
  I_a :=0;
  for i in 1:nLfun-1 loop
    if datafun[i].z < 0 then
      I_a :=I_a + 0.5*mol_i[i]*datafun[i].z^2;
    end if;
  end for;
end calc_I_a;
