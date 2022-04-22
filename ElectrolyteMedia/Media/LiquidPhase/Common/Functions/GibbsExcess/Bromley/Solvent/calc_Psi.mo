within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_Psi "calculates Psi from Zemaitis1986,page 238"

  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] Psi;
protected
  Real[nLfun - 1,nLfun - 1] aI=calc_aI(X);
algorithm

  for i in 1:nLfun-1 loop
    for j in 1:nLfun-1 loop
      if datafun[i].z > 0 and datafun[j].z < 0 then
        Psi[i, j] := 2/(aI[i, j] + 1e-25)*((1 + 2*aI[i, j])/(1 + aI[i, j])^2 -
          log(1 + aI[i, j])/(aI[i, j] + 1e-25));
      end if;
    end for;
  end for;

end calc_Psi;
