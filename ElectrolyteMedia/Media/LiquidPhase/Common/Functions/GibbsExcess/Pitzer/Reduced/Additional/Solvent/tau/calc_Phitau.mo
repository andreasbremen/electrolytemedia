within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solvent.tau;
function calc_Phitau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real Phi;

protected
  Real m_S;
  Real A_Phi_lntau = GibbsExcess.DebyeHueckel.Reduced.Additional.dg_dtau.calc_dAln_dtau(T,p);
  Real I = calc_I(X);
  Real[nLifun] mol_i = calc_mfromX(X);

  Real Z = GibbsExcess.Pitzer.Solute.calc_Z(X);
  Real[nLifun,nLifun] C_ijtau = GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_C_ijtau(T);
  Real[nLifun,nLifun] B_Phi_ijtau = calc_B_Phi_ijtau(T,X);
  Real[nLifun,nLifun] Phi_Phi_ijtau = calc_Phi_Phi_ijtau(T,p,X);
  Real [nLifun,nLifun,nLifun] psi_ijktau = GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_psi_ijktau(T);
  Real [nLifun,nLifun] lambda_ijtau = GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_lambda_ijtau(T);
  Real [nLifun,nLifun,nLifun] zeta_ijktau = GibbsExcess.Pitzer.Reduced.Additional.Solute.tau.calc_zeta_ijktau(T);
algorithm
//   for i in 1:nLifun loop
//     if datafun[i].z <> 0 then
//       m_S :=m_S + mol_i[i];
//     end if;
//   end for;

  m_S :=sum(mol_i);

  Phi := -(A_Phi_lntau*I^(3/2))/(1 + 1.2*I^(1/2));

  //cation-anion
  for i in 1:nLifun loop
    if datafun[i].z > 0 then
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          Phi :=Phi + mol_i[i]*mol_i[j]*(B_Phi_ijtau[i,j] + Z*C_ijtau[i,j]);
        end if;
      end for;
    end if;
  end for;

  //cation-cation
  for i in 1:nLifun loop
    if datafun[i].z > 0 then
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          Phi :=Phi + mol_i[i]*mol_i[j]*(Phi_Phi_ijtau[i,j]);
          //cation-cation-anion
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              Phi :=Phi + mol_i[i]*mol_i[j]*mol_i[k]*psi_ijktau[i, j, k];
            end if;
          end for;
        end if;
      end for;
    end if;
  end for;

  //anion-anion
  for i in 1:nLifun loop
    if datafun[i].z < 0 then
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          Phi :=Phi + mol_i[i]*mol_i[j]*(Phi_Phi_ijtau[i,j]);
          //anion-anion-cation
          for k in 1:nLifun loop
            if datafun[k].z > 0 then
              Phi :=Phi + mol_i[i]*mol_i[j]*mol_i[k]*psi_ijktau[i, j, k];
            end if;
          end for;
        end if;
      end for;
    end if;
  end for;

  //neutral-neutral
  for i in 1:nLifun loop
    if datafun[i].z == 0 then
      Phi :=Phi + 0.5*mol_i[i]^2*lambda_ijtau[i, i];
      for j in 1:nLifun loop
        if datafun[j].z == 0 and i <> j then
          Phi :=Phi + mol_i[i]*mol_i[j]*lambda_ijtau[i, j];
        end if;
      end for;
    end if;
  end for;

  //charged-neutral
  for i in 1:nLifun loop
    if datafun[i].z == 0 then
      for j in 1:nLifun loop
        if datafun[j].z <> 0 then
          Phi :=Phi + mol_i[i]*mol_i[j]*lambda_ijtau[i, j];
        end if;
        if datafun[j].z > 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              Phi :=Phi + mol_i[i]*mol_i[j]*mol_i[k]*zeta_ijktau[i, j, k];
            end if;
          end for;
        end if;
      end for;
    end if;
  end for;

  Phi := 2/m_S*Phi;

end calc_Phitau;
