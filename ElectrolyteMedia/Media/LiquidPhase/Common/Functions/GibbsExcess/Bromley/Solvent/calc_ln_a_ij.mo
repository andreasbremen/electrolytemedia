within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_ln_a_ij
  "calculates natural logarithm of water activity from Bromley1973"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] ln_a_s;
protected
  Real mol_i[nLfun-1] = calc_mfromX(X);
  parameter SI.MolarMass M_s = 0.0180152;
  Real[nLfun - 1,nLfun - 1] Phi_i=calc_Phi_ij(T,p,X);
algorithm
  for i in 1:nLfun-1 loop
    for j in 1:nLfun-1 loop
      if datafun[i].z > 0 and datafun[j].z < 0 then
        ln_a_s[i,j] := -IF97.MH2O *Phi_i[i,j]*(mol_i[i] + mol_i[j]);
      end if;
    end for;
  end for;
end calc_ln_a_ij;
