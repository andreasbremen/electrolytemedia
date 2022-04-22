within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tautau;
function calc_loggammatautau "Helper function to calculate Gibbs derivative"

  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun] log_gammatautau;

protected
  Real[nLifun] m_i=calc_mfromX(X);
  Real[nLifun,nLifun] lambda_ijtautau=calc_lambda_ijtautau(T);
  Real[nLifun,nLifun,nLifun] psi_ijktautau=calc_psi_ijktautau(T);
  Real[nLifun,nLifun,nLifun] zeta_ijktautau=calc_zeta_ijktautau(T);
  Real Ftautau=calc_Ftautau(
      T,
      p,
      X);
  Real[nLifun,nLifun] C_ijtautau=calc_C_ijtautau(T);
  Real[nLifun,nLifun] B_ijtautau=calc_B_ijtautau(T, X);
  Real[nLifun,nLifun] Phi_ijtautau=calc_Phi_ijtautau(
      T,
      p,
      X);
  Real Z=Pitzer.Solute.calc_Z(X);
  Real[nLifun] ln_gammatautau;
  Real C_ijtautau_temp;
algorithm
  for i in 1:nLifun loop
    C_ijtautau_temp :=0;
    //cations
    if datafun[i].z > 0 then
      //solvent contribution
      ln_gammatautau[i] :=datafun[i].z^2*Ftautau;
      //counter-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          ln_gammatautau[i] :=ln_gammatautau[i] + m_i[j]*(2*B_ijtautau[i, j] + Z*C_ijtautau[i, j]);
        end if;
      end for;
      for j in 1:nLifun loop
        for k in j:nLifun loop
          C_ijtautau_temp :=C_ijtautau_temp + m_i[j]*m_i[k]*C_ijtautau[j, k];
        end for;
      end for;
      ln_gammatautau[i] :=ln_gammatautau[i] + abs(datafun[i].z)*C_ijtautau_temp;
      //like-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z > 0  and j <> i then
          ln_gammatautau[i] :=ln_gammatautau[i] + m_i[j]*2*Phi_ijtautau[i, j];
        end if;
      end for;
      //2 cations, 1 anion
      for j in 1:nLifun loop
        if datafun[j].z > 0  and j <> i then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              ln_gammatautau[i] :=ln_gammatautau[i] + m_i[j]*m_i[k]*psi_ijktautau[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //2 anions, 1 cation
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0  and j <> k then
              ln_gammatautau[i] :=ln_gammatautau[i] + m_i[j]*m_i[k]*psi_ijktautau[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //cation-neutral
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          ln_gammatautau[i] :=ln_gammatautau[i] + 2*m_i[j]*lambda_ijtautau[i, j];
        end if;
      end for;
      //1 cation, 1 neutral, 1 anion
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              ln_gammatautau[i] :=ln_gammatautau[i] + 6*m_i[j]*m_i[k]*zeta_ijktautau[i, j, k];
            end if;
          end for;
        end if;
      end for;

    //anions
    elseif datafun[i].z < 0 then
      //solvent contribution
      ln_gammatautau[i] :=datafun[i].z^2*Ftautau;
      //counter-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          ln_gammatautau[i] :=ln_gammatautau[i] + m_i[j]*(2*B_ijtautau[i, j] + Z*C_ijtautau[i, j]);
        end if;
      end for;
      for j in 1:nLifun loop
        for k in j:nLifun loop
          C_ijtautau_temp :=C_ijtautau_temp + m_i[j]*m_i[k]*C_ijtautau[j, k];
        end for;
      end for;
      ln_gammatautau[i] :=ln_gammatautau[i] + abs(datafun[i].z)*C_ijtautau_temp;
      //like-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z < 0  and j <> i then
          ln_gammatautau[i] :=ln_gammatautau[i] + m_i[j]*2*Phi_ijtautau[i, j];
        end if;
      end for;
      //2 anions, 1 cation
      for j in 1:nLifun loop
        if datafun[j].z < 0  and j <> i then
          for k in 1:nLifun loop
            if datafun[k].z > 0 then
              ln_gammatautau[i] :=ln_gammatautau[i] + m_i[j]*m_i[k]*psi_ijktautau[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //2 cations, 1 anion
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          for k in 1:nLifun loop
            if datafun[k].z > 0  and j <> k then
              ln_gammatautau[i] :=ln_gammatautau[i] + m_i[j]*m_i[k]*psi_ijktautau[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //anion-neutral
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          ln_gammatautau[i] :=ln_gammatautau[i] + 2*m_i[j]*lambda_ijtautau[i, j];
        end if;
      end for;
      //1 anion, 1 neutral, 1 cation
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          for k in 1:nLifun loop
            if datafun[k].z > 0 then
              ln_gammatautau[i] :=ln_gammatautau[i] + 6*m_i[j]*m_i[k]*zeta_ijktautau[i, j, k];
            end if;
          end for;
        end if;
      end for;

    //neutrals
    else
      //neutral-neutral
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          ln_gammatautau[i] := ln_gammatautau[i] + 2*m_i[j]*lambda_ijtautau[i,j];
        end if;
      end for;
      //neutral-cation
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          ln_gammatautau[i] :=ln_gammatautau[i] + 2*m_i[j]*lambda_ijtautau[i, j];
        end if;
      end for;
      //neutral-anion
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          ln_gammatautau[i] :=ln_gammatautau[i] + 2*m_i[j]*lambda_ijtautau[i, j];
        end if;
      end for;
      //neutral-cation-anion
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              ln_gammatautau[i] :=ln_gammatautau[i] + m_i[j]*m_i[k]*zeta_ijktautau[i, j, k];
            end if;
          end for;
        end if;
      end for;

    end if;
  end for;

  log_gammatautau :=ln_gammatautau/log(10);
end calc_loggammatautau;
