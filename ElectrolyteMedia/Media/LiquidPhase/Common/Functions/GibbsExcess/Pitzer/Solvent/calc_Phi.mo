within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solvent;
function calc_Phi
  "Calculates osmotic coefficient based on Pitzer model"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real Phi;

protected
  Real m_S;
  Real A_Phi_ln = DebyeHueckel.calc_A_ln(T,p);
  Real I = calc_I(X);
  Real[nLifun] mol_i = calc_mfromX(X);

  Real Z = Solute.calc_Z(X);
  Real[nLifun,nLifun] C_ij = Solute.calc_C_ij(T);
  Real[nLifun,nLifun] B_Phi_ij = calc_B_Phi_ij(T,X);
  Real[nLifun,nLifun] Phi_Phi_ij = calc_Phi_Phi_ij(T,p,X);
  Real [nLifun,nLifun,nLifun] psi_ijk = Solute.calc_psi_ijk(T);
  Real [nLifun,nLifun] lambda_ij = Solute.calc_lambda_ij(T);
  Real [nLifun,nLifun,nLifun] zeta_ijk = Solute.calc_zeta_ijk(T);
algorithm
//   for i in 1:nLifun loop
//     if datafun[i].z <> 0 then
//       m_S :=m_S + mol_i[i];
//     end if;
//   end for;

  m_S :=sum(mol_i);

  Phi := -(A_Phi_ln*I^(3/2))/(1 + 1.2*I^(1/2));

  //cation-anion
  for i in 1:nLifun loop
    if datafun[i].z > 0 then
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          Phi :=Phi + mol_i[i]*mol_i[j]*(B_Phi_ij[i,j] + Z*C_ij[i,j]);
        end if;
      end for;
    end if;
  end for;

  //cation-cation
  for i in 1:nLifun loop
    if datafun[i].z > 0 then
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          Phi :=Phi + mol_i[i]*mol_i[j]*(Phi_Phi_ij[i,j]);
          //cation-cation-anion
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              Phi :=Phi + mol_i[i]*mol_i[j]*mol_i[k]*psi_ijk[i, j, k];
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
          Phi :=Phi + mol_i[i]*mol_i[j]*(Phi_Phi_ij[i,j]);
          //anion-anion-cation
          for k in 1:nLifun loop
            if datafun[k].z > 0 then
              Phi :=Phi + mol_i[i]*mol_i[j]*mol_i[k]*psi_ijk[i, j, k];
            end if;
          end for;
        end if;
      end for;
    end if;
  end for;

  //neutral-neutral
  for i in 1:nLifun loop
    if datafun[i].z == 0 then
      Phi :=Phi + 0.5*mol_i[i]^2*lambda_ij[i, i];
      for j in 1:nLifun loop
        if datafun[j].z == 0 and i <> j then
          Phi :=Phi + mol_i[i]*mol_i[j]*lambda_ij[i, j];
        end if;
      end for;
    end if;
  end for;

  //charged-neutral
  for i in 1:nLifun loop
    if datafun[i].z == 0 then
      for j in 1:nLifun loop
        if datafun[j].z <> 0 then
          Phi :=Phi + mol_i[i]*mol_i[j]*lambda_ij[i, j];
        end if;
        if datafun[j].z > 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              Phi :=Phi + mol_i[i]*mol_i[j]*mol_i[k]*zeta_ijk[i, j, k];
            end if;
          end for;
        end if;
      end for;
    end if;
  end for;

  if m_S > 0 then
    Phi :=1 + 2/m_S*Phi;
  else
    Phi :=1;
  end if;

end calc_Phi;
