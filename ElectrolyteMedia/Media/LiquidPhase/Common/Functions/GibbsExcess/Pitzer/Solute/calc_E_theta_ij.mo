within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_E_theta_ij "Unsymmetrical mixing, Pitzer p.123"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] E_theta_ij;
protected
  Real[nLifun,nLifun] x_ij=Solute.calc_x_ij(
      T,
      p,
      X);
  Real I = calc_I(X);
  Real J_ij;
  Real J_ii;
  Real J_jj;
algorithm
  if I > 0 then
    for i in 1:nLifun loop
      for j in 1:nLifun loop
        if i <> j then
          if datafun[i].z * datafun[j].z > 0 and datafun[i].z <> datafun[j].z then
            J_ij := Solute.calc_J(x_ij[i, j]);
            J_ii := Solute.calc_J(x_ij[i, i]);
            J_jj := Solute.calc_J(x_ij[j, j]);
            E_theta_ij[i,j] :=datafun[i].z*datafun[j].z/(4*I)*(J_ij - J_ii/2 - J_jj/2);
          end if;
        end if;
      end for;
    end for;
  end if;

end calc_E_theta_ij;
