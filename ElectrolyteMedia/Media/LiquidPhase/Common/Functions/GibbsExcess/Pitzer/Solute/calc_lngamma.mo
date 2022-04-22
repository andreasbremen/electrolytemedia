within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_lngamma
  "Calculates natural logarithm of activity coefficient in molality base of aqueous species from Pitzer model"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun] ln_gamma;
protected
  Real [nLifun] m_i = calc_mfromX(X);
  Real[nLifun,nLifun] beta0=Solute.calc_beta0_ij(T);
  Real[nLifun,nLifun] beta1=Solute.calc_beta1_ij(T);
  Real[nLifun,nLifun] beta2=Solute.calc_beta2_ij(T);
  Real[nLifun,nLifun] c_ij=Solute.calc_c_ij(T);
  Real[nLifun,nLifun] theta_ij=Solute.calc_theta_ij(T);
  Real[nLifun,nLifun] lambda_ij=Solute.calc_lambda_ij(T);
  Real[nLifun,nLifun,nLifun] psi_ijk=Solute.calc_psi_ijk(T);
  Real[nLifun,nLifun,nLifun] zeta_ijk=Solute.calc_zeta_ijk(T);
  Real F=Solute.calc_F(
      T,
      p,
      X);
  Real[nLifun,nLifun] C_ij=Solute.calc_C_ij(T);
  Real[nLifun,nLifun] B_ij=Solute.calc_B_ij(T, X);
  Real[nLifun,nLifun] Phi_ij=Solute.calc_Phi_ij(
      T,
      p,
      X);
  Real Z=Solute.calc_Z(X);
  Real C_ij_temp;
algorithm
  for i in 1:nLifun loop
    C_ij_temp :=0;
    //cations
    if datafun[i].z > 0 then
      //solvent contribution
      ln_gamma[i] :=datafun[i].z^2*F;
      //counter-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          ln_gamma[i] :=ln_gamma[i] + m_i[j]*(2*B_ij[i, j] + Z*C_ij[i, j]);
        end if;
      end for;
      for j in 1:nLifun loop
        for k in j:nLifun loop
          C_ij_temp :=C_ij_temp + m_i[j]*m_i[k]*C_ij[j, k];
        end for;
      end for;
      ln_gamma[i] :=ln_gamma[i] + abs(datafun[i].z)*C_ij_temp;
      //like-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z > 0  and j <> i then
          ln_gamma[i] :=ln_gamma[i] + m_i[j]*2*Phi_ij[i, j];
        end if;
      end for;
      //2 cations, 1 anion
      for j in 1:nLifun loop
        if datafun[j].z > 0  and j <> i then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              ln_gamma[i] :=ln_gamma[i] + m_i[j]*m_i[k]*psi_ijk[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //2 anions, 1 cation
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0  and j <> k then
              ln_gamma[i] :=ln_gamma[i] + m_i[j]*m_i[k]*psi_ijk[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //cation-neutral
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          ln_gamma[i] :=ln_gamma[i] + 2*m_i[j]*lambda_ij[i, j];
        end if;
      end for;
      //1 cation, 1 neutral, 1 anion
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              ln_gamma[i] :=ln_gamma[i] + 6*m_i[j]*m_i[k]*zeta_ijk[i, j, k];
            end if;
          end for;
        end if;
      end for;

    //anions
    elseif datafun[i].z < 0 then
      //solvent contribution
      ln_gamma[i] :=datafun[i].z^2*F;
      //counter-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          ln_gamma[i] :=ln_gamma[i] + m_i[j]*(2*B_ij[i, j] + Z*C_ij[i, j]);
        end if;
      end for;
      for j in 1:nLifun loop
        for k in j:nLifun loop
          C_ij_temp :=C_ij_temp + m_i[j]*m_i[k]*C_ij[j, k];
        end for;
      end for;
      ln_gamma[i] :=ln_gamma[i] + abs(datafun[i].z)*C_ij_temp;
      //like-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z < 0  and j <> i then
          ln_gamma[i] :=ln_gamma[i] + m_i[j]*2*Phi_ij[i, j];
        end if;
      end for;
      //2 anions, 1 cation
      for j in 1:nLifun loop
        if datafun[j].z < 0  and j <> i then
          for k in 1:nLifun loop
            if datafun[k].z > 0 then
              ln_gamma[i] :=ln_gamma[i] + m_i[j]*m_i[k]*psi_ijk[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //2 cations, 1 anion
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          for k in 1:nLifun loop
            if datafun[k].z > 0  and j <> k then
              ln_gamma[i] :=ln_gamma[i] + m_i[j]*m_i[k]*psi_ijk[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //anion-neutral
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          ln_gamma[i] :=ln_gamma[i] + 2*m_i[j]*lambda_ij[i, j];
        end if;
      end for;
      //1 anion, 1 neutral, 1 cation
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          for k in 1:nLifun loop
            if datafun[k].z > 0 then
              ln_gamma[i] :=ln_gamma[i] + 6*m_i[j]*m_i[k]*zeta_ijk[i, j, k];
            end if;
          end for;
        end if;
      end for;

    //neutrals
    else
      //neutral-neutral
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          ln_gamma[i] := ln_gamma[i] + 2*m_i[j]*lambda_ij[i,j];
        end if;
      end for;
      //neutral-cation
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          ln_gamma[i] :=ln_gamma[i] + 2*m_i[j]*lambda_ij[i, j];
        end if;
      end for;
      //neutral-anion
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          ln_gamma[i] :=ln_gamma[i] + 2*m_i[j]*lambda_ij[i, j];
        end if;
      end for;
      //neutral-cation-anion
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              ln_gamma[i] :=ln_gamma[i] + m_i[j]*m_i[k]*zeta_ijk[i, j, k];
            end if;
          end for;
        end if;
      end for;

    end if;
  end for;

 annotation(smoothOrder(zeroDerivative=data)=2);
end calc_lngamma;
