within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_C_ij "Calculates cation i - anion j interaction C_ij"
  input Modelica.SIunits.Temperature T;
  output Real [nLifun, nLifun] C_ij;
protected
  Real [nLifun,nLifun] c_ij = calc_c_ij(T);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if abs(datafun[i].z) > 0 and abs(datafun[j].z) > 0 then
        C_ij [i,j] := c_ij[i,j]/(2*(abs(datafun[i].z*datafun[j].z))^0.5);
      end if;
    end for;
  end for;
end calc_C_ij;

function calc_c_ij
  "Calculates temperature dependent Pitzer polynomial for binary cation-anion interaction"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] c_ij;
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      c_ij[i,j] := interactionfun.c_ij[i,j,1] + interactionfun.c_ij[i,j,2]*(1/T-1/Tref) + interactionfun.c_ij[i,j,3]*log(T/Tref) + interactionfun.c_ij[i,j,4]*(T-Tref) + interactionfun.c_ij[i,j,5]*(T^2-Tref^2) + interactionfun.c_ij[i,j,6]*(1/T^2-1/Tref^2) + interactionfun.c_ij[i,j,7]*T^3 + interactionfun.c_ij[i,j,8]*1/(T-263);
    end for;
  end for;
//   c_ij := interactionfun.c_ij[:,:,1] + interactionfun.c_ij[:,:,2]*(1/T-1/Tref) + interactionfun.c_ij[:,:,3]*log(T/Tref) + interactionfun.c_ij[:,:,4]*(T-Tref) + interactionfun.c_ij[:,:,5]*(T^2-Tref^2) + interactionfun.c_ij[:,:,6]*(1/T^2-1/Tref^2) + interactionfun.c_ij[:,:,7]*T^3 + interactionfun.c_ij[:,:,8]*1/(T-263);

end calc_c_ij;
