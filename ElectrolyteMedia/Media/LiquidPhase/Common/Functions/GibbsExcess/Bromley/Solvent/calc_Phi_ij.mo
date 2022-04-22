within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Solvent;
function calc_Phi_ij "calculate osmotic coefficient for one electrolyte"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] Phi_i;
protected
  Real I = calc_I(X);
  Real A=DebyeHueckel.calc_A_log(T,p);
  Real sigma=calc_sigma(X);
  Real[nLfun - 1,nLfun - 1] Psi=calc_Psi(X);
algorithm

  for i in 1:nLfun-1 loop
    for j in 1:nLfun-1 loop
      if datafun[i].z > 0 and datafun[j].z < 0 then
        Phi_i[i,j] := 1 - log(10)*(A*abs(datafun[i].z*datafun[j].z)*I^0.5/3*sigma - (0.06 +
          0.6*interactionfun.Bromley_ij[i, j])*abs(datafun[i].z*datafun[j].z)*I/2*Psi[i,j] -
          interactionfun.Bromley_ij[i, j]*I/2);
      end if;
    end for;
  end for;
end calc_Phi_ij;
