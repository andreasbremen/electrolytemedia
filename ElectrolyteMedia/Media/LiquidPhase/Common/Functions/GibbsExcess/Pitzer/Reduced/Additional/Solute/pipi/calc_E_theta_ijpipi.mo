within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pipi;
function calc_E_theta_ijpipi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] E_theta_ijpipi;
protected
  Real[nLifun,nLifun] x_ij=GibbsExcess.Pitzer.Solute.calc_x_ij(T,p,X);
  Real[nLifun,nLifun] x_ijpi = pi.calc_x_ijpi(T,p,X);
  Real[nLifun,nLifun] x_ijpipi = calc_x_ijpipi(T,p,X);
  Real I = calc_I(X);
  Real J_ijdx;
  Real J_iidx;
  Real J_jjdx;
  Real J_ijdxdx;
  Real J_iidxdx;
  Real J_jjdxdx;

algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if i <> j then
        if datafun[i].z*datafun[j].z > 0 and datafun[i].z <> datafun[j].z then
          J_ijdx :=GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[i, j]);
          J_iidx :=GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[i, i]);
          J_jjdx :=GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[j, j]);
          J_ijdxdx :=tau.calc_J_dxdx(x_ij[i, j]);
          J_iidxdx :=tau.calc_J_dxdx(x_ij[i, i]);
          J_jjdxdx :=tau.calc_J_dxdx(x_ij[j, j]);
          E_theta_ijpipi[i,j] :=datafun[i].z*datafun[j].z/(4*I)*(J_ijdxdx*x_ijpi[i,j]^2+J_ijdx*x_ijpipi[i,j]-1/2*(J_iidxdx*x_ijpi[i,i]^2+J_iidx*x_ijpipi[i,i])-1/2*(J_jjdxdx*x_ijpi[j,j]^2+J_jjdx*x_ijpipi[j,j]));
        end if;
      end if;
    end for;
  end for;

end calc_E_theta_ijpipi;
