within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.pi;
function calc_E_theta_ijpi "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun,nLifun] E_theta_ijpi;
protected
  Real[nLifun,nLifun] x_ij=GibbsExcess.Pitzer.Solute.calc_x_ij(T,p,X);
  Real[nLifun,nLifun] x_ijpi = calc_x_ijpi(T,p,X);
  Real I = calc_I(X);
  Real J_ijdx;
  Real J_iidx;
  Real J_jjdx;

algorithm
  for i in 1:nLifun loop
    for j in 1:nLifun loop
      if i <> j then
        if datafun[i].z*datafun[j].z > 0 and datafun[i].z <> datafun[j].z then
          J_ijdx :=GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[i, j]);
          J_iidx :=GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[i, i]);
          J_jjdx :=GibbsExcess.Pitzer.Solute.calc_J_dx(x_ij[j, j]);
          E_theta_ijpi[i,j] :=datafun[i].z*datafun[j].z/(4*I)*(J_ijdx*x_ijpi[i,j]-1/2*J_iidx*x_ijpi[i,i]-1/2*J_jjdx*x_ijpi[j,j]);
        end if;
      end if;
    end for;
  end for;

end calc_E_theta_ijpi;
