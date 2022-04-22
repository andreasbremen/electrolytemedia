within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_der_E_theta_ij "Unsymmetrical mixing, Pitzer p.124"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] der_E_theta_ij;
protected
  Real[nLifun,nLifun] x_ij=Solute.calc_x_ij(
      T,
      p,
      X);
  Real[nLifun,nLifun] E_theta_ij=Solute.calc_E_theta_ij(
      T,
      p,
      X);
  Real I = calc_I(X);
  Real J_ij_dx;
  Real J_ii_dx;
  Real J_jj_dx;
algorithm
  if I > 0 then
    for i in 1:nLifun loop
      for j in 1:nLifun loop
        if i <> j then
          if datafun[i].z * datafun[j].z > 0 and datafun[i].z <> datafun[j].z then
            J_ij_dx := Solute.calc_J_dx(x_ij[i, j]);
            J_ii_dx := Solute.calc_J_dx(x_ij[i, i]);
            J_jj_dx := Solute.calc_J_dx(x_ij[j, j]);
            der_E_theta_ij[i,j] :=-E_theta_ij[i,j]/I + datafun[i].z*datafun[j].z/(8*I^2)*(x_ij[i,j]*J_ij_dx - x_ij[i,i]*J_ii_dx/2 - x_ij[j,j]*J_jj_dx/2);
          end if;
        end if;
      end for;
    end for;
  end if;
end calc_der_E_theta_ij;
