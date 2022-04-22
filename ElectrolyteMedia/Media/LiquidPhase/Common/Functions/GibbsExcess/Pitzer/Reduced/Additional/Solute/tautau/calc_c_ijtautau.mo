within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_c_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  output Real[nLifun,nLifun] c_ij;
protected
  Real tau = calc_tau(T);
algorithm
  c_ij := 2*Tref/tau^3*(-interactionfun.c_ij[:,:,2]*(1/T^2) + interactionfun.c_ij[:,:,3]*1/T + interactionfun.c_ij[:,:,4] + interactionfun.c_ij[:,:,5]*(2*T) + interactionfun.c_ij[:,:,6]*(-2/T^3) + interactionfun.c_ij[:,:,7]*3*T^2 + interactionfun.c_ij[:,:,8]*(-1)/(T-263)^(2))
              + (-Tref/tau^2)^2*(2*interactionfun.c_ij[:,:,2]*(1/T^3) - interactionfun.c_ij[:,:,3]*1/T^2 + interactionfun.c_ij[:,:,5]*2 + interactionfun.c_ij[:,:,6]*(6/T^4) + interactionfun.c_ij[:,:,7]*6*T + 2*interactionfun.c_ij[:,:,8]*1/(T-263)^(3));

end calc_c_ijtautau;

function calc_C_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  output Real [nLifun, nLifun] C_ijtautau;
protected
  Real[nLifun,nLifun] c_ijtautau=calc_c_ijtautau(
                                           T);
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if abs(datafun[i].z) > 0 and abs(datafun[j].z) > 0 then
        C_ijtautau [i,j] := c_ijtautau[i,j]/(2*(abs(datafun[i].z*datafun[j].z))^0.5);
      end if;
    end for;
  end for;
end calc_C_ijtautau;
