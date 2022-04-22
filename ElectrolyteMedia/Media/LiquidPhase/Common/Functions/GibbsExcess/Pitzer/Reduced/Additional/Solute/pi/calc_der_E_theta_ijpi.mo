within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pi;
function calc_der_E_theta_ijpi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] der_E_theta_ijpi;
protected
  Real[nLifun,nLifun] x_ij=GibbsExcess.Pitzer.Solute.calc_x_ij(
      T,
      p,
      X);
  Real[nLifun,nLifun] x_ijpi = calc_x_ijpi(T,p,X);
  Real[nLifun,nLifun] E_theta_ijpi=calc_E_theta_ijpi(
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
algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if i <> j then
        if datafun[i].z*datafun[j].z > 0 and datafun[i].z <> datafun[j].z then
          J_ij_dx := GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[i, j]);
          J_ii_dx := GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[i, i]);
          J_jj_dx := GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[j, j]);
          J_ij_dxdx := GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_J_dxdx(x_ij[i, j]);
          J_ii_dxdx := GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_J_dxdx(x_ij[i, i]);
          J_jj_dxdx := GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_J_dxdx(x_ij[j, j]);
          der_E_theta_ijpi[i,j] :=-E_theta_ijpi[i,j]/I + datafun[i].z*datafun[j].z/(8*I^2)*(x_ijpi[i,j]*J_ij_dx+x_ij[i,j]*J_ij_dxdx*x_ijpi[i,j]-1/2*(x_ijpi[i,i]*J_ii_dx+x_ij[i,i]*J_ii_dxdx*x_ijpi[i,i])-1/2*(x_ijpi[j,j]*J_jj_dx+x_ij[j,j]*J_jj_dxdx*x_ijpi[j,j]));
        end if;
      end if;
    end for;
  end for;

end calc_der_E_theta_ijpi;
