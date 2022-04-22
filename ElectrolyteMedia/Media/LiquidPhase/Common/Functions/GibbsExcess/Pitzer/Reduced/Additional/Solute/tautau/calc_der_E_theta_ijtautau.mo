within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_der_E_theta_ijtautau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] der_E_theta_ijtautau;
protected
  Real[nLifun,nLifun] x_ij=Pitzer.Solute.calc_x_ij(
      T,
      p,
      X);
  Real[nLifun,nLifun] x_ijtau=tau.calc_x_ijtau(
      T,
      p,
      X);
      Real[nLifun,nLifun] x_ijtautau=calc_x_ijtautau(
      T,
      p,
      X);
  Real[nLifun,nLifun] E_theta_ijtautau=calc_E_theta_ijtautau(
      T,
      p,
      X);
  Real I = calc_I(X);
  Real J_ij_dx;
  Real J_ii_dx;
  Real J_jj_dx;
  Real J_ij_dxdx;
  Real J_ii_dxdx;
  Real J_jj_dxdx;
  Real J_ij_dxdxdx;
  Real J_ii_dxdxdx;
  Real J_jj_dxdxdx;
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if i <> j then
        if datafun[i].z * datafun[j].z > 0 and datafun[i].z <> datafun[j].z then
          J_ij_dx := Pitzer.Solute.calc_J_dx(x_ij[i, j]);
          J_ii_dx := Pitzer.Solute.calc_J_dx(x_ij[i, i]);
          J_jj_dx := Pitzer.Solute.calc_J_dx(x_ij[j, j]);
          J_ij_dxdx := tau.calc_J_dxdx(x_ij[i, j]);
          J_ii_dxdx := tau.calc_J_dxdx(x_ij[i, i]);
          J_jj_dxdx := tau.calc_J_dxdx(x_ij[j, j]);
          J_ij_dxdxdx := calc_J_dxdxdx(x_ij[i, j]);
          J_ii_dxdxdx := calc_J_dxdxdx(x_ij[i, i]);
          J_jj_dxdxdx := calc_J_dxdxdx(x_ij[j, j]);
          der_E_theta_ijtautau[i,j] :=-E_theta_ijtautau[i,j]/I + datafun[i].z*datafun[j].z/(8*I^2)*(x_ijtautau[i,j]*J_ij_dx+2*x_ijtau[i,j]^2*J_ij_dxdx+x_ij[i,j]*J_ij_dxdxdx*x_ijtau[i,j]^2+x_ij[i,j]*J_ij_dxdx*x_ijtautau[i,j]-1/2*((x_ijtautau[i,i]*J_ii_dx+2*x_ijtau[i,i]^2*J_ii_dxdx+x_ij[i,i]*J_ii_dxdxdx*x_ijtau[i,i]^2+x_ij[i,i]*J_ii_dxdx*x_ijtautau[i,i]))-1/2*((x_ijtautau[j,j]*J_jj_dx+2*x_ijtau[j,j]^2*J_jj_dxdx+x_ij[j,j]*J_jj_dxdxdx*x_ijtau[j,j]^2+x_ij[j,j]*J_jj_dxdx*x_ijtautau[j,j])));
        end if;
      end if;
    end for;
  end for;

end calc_der_E_theta_ijtautau;
