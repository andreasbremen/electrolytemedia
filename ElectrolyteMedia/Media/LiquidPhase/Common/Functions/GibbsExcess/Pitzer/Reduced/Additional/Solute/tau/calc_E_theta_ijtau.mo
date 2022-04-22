within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_E_theta_ijtau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] E_theta_ijtau;
protected
  Real[nLifun,nLifun] x_ij=Pitzer.Solute.calc_x_ij(T,p,X);
  Real[nLifun,nLifun] x_ijtau=calc_x_ijtau(T,p,X);
  Real I = calc_I(X);
  Real J_ij_dx;
  Real J_ii_dx;
  Real J_jj_dx;

algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if i <> j then
        if datafun[i].z*datafun[j].z > 0 and datafun[i].z <> datafun[j].z then
          J_ij_dx := GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[i, j]);
          J_ii_dx := GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[i, i]);
          J_jj_dx := GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[j, j]);
          E_theta_ijtau[i,j] :=datafun[i].z*datafun[j].z/(4*I)*(J_ij_dx *x_ijtau[i,j] - J_ii_dx/2 *x_ijtau[i,i] - J_jj_dx/2 *x_ijtau[j,j]);
        end if;
      end if;
    end for;
  end for;

end calc_E_theta_ijtau;
