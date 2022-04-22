within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Reduced.Additional.Solute.tau;
function calc_loggammatau "Helper function to calculate Gibbs derivative"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real[nLifun] log_gammatau;

protected
  Real [nLifun] m_i = calc_mfromX(X);
  Real[nLifun,nLifun] lambda_ijtau=calc_lambda_ijtau(T);
  Real[nLifun,nLifun,nLifun] psi_ijktau=calc_psi_ijktau(T);
  Real[nLifun,nLifun,nLifun] zeta_ijktau=calc_zeta_ijktau(T);
  Real Ftau=calc_Ftau(T,p,X);
  Real[nLifun,nLifun] C_ijtau=calc_C_ijtau(T);
  Real[nLifun,nLifun] B_ijtau=calc_B_ijtau(T, X);
  Real[nLifun,nLifun] Phi_ijtau=calc_Phi_ijtau(T,p,X);
  Real Z=Pitzer.Solute.calc_Z(X);
  Real[nLifun] ln_gammatau;
  Real C_ijtau_temp;
algorithm
  for i in 1:nLifun loop
    C_ijtau_temp :=0;
    //cations
    if datafun[i].z > 0 then
      //solvent contribution
      ln_gammatau[i] :=datafun[i].z^2*Ftau;
      //counter-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          ln_gammatau[i] :=ln_gammatau[i] + m_i[j]*(2*B_ijtau[i, j] + Z*C_ijtau[i, j]);
        end if;
      end for;
      for j in 1:nLifun loop
        for k in j:nLifun loop
          C_ijtau_temp :=C_ijtau_temp + m_i[j]*m_i[k]*C_ijtau[j, k];
        end for;
      end for;
      ln_gammatau[i] :=ln_gammatau[i] + abs(datafun[i].z)*C_ijtau_temp;
      //like-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z > 0  and j <> i then
          ln_gammatau[i] :=ln_gammatau[i] + m_i[j]*2*Phi_ijtau[i, j];
        end if;
      end for;
      //2 cations, 1 anion
      for j in 1:nLifun loop
        if datafun[j].z > 0  and j <> i then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              ln_gammatau[i] :=ln_gammatau[i] + m_i[j]*m_i[k]*psi_ijktau[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //2 anions, 1 cation
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0  and j <> k then
              ln_gammatau[i] :=ln_gammatau[i] + m_i[j]*m_i[k]*psi_ijktau[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //cation-neutral
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          ln_gammatau[i] :=ln_gammatau[i] + 2*m_i[j]*lambda_ijtau[i, j];
        end if;
      end for;
      //1 cation, 1 neutral, 1 anion
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              ln_gammatau[i] :=ln_gammatau[i] + 6*m_i[j]*m_i[k]*zeta_ijktau[i, j, k];
            end if;
          end for;
        end if;
      end for;

    //anions
    elseif datafun[i].z < 0 then
      //solvent contribution
      ln_gammatau[i] :=datafun[i].z^2*Ftau;
      //counter-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          ln_gammatau[i] :=ln_gammatau[i] + m_i[j]*(2*B_ijtau[i, j] + Z*C_ijtau[i, j]);
        end if;
      end for;
      for j in 1:nLifun loop
        for k in j:nLifun loop
          C_ijtau_temp :=C_ijtau_temp + m_i[j]*m_i[k]*C_ijtau[j, k];
        end for;
      end for;
      ln_gammatau[i] :=ln_gammatau[i] + abs(datafun[i].z)*C_ijtau_temp;
      //like-charged contribution
      for j in 1:nLifun loop
        if datafun[j].z < 0  and j <> i then
          ln_gammatau[i] :=ln_gammatau[i] + m_i[j]*2*Phi_ijtau[i, j];
        end if;
      end for;
      //2 anions, 1 cation
      for j in 1:nLifun loop
        if datafun[j].z < 0  and j <> i then
          for k in 1:nLifun loop
            if datafun[k].z > 0 then
              ln_gammatau[i] :=ln_gammatau[i] + m_i[j]*m_i[k]*psi_ijktau[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //2 cations, 1 anion
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          for k in 1:nLifun loop
            if datafun[k].z > 0  and j <> k then
              ln_gammatau[i] :=ln_gammatau[i] + m_i[j]*m_i[k]*psi_ijktau[i, j, k];
            end if;
          end for;
        end if;
      end for;
      //anion-neutral
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          ln_gammatau[i] :=ln_gammatau[i] + 2*m_i[j]*lambda_ijtau[i, j];
        end if;
      end for;
      //1 anion, 1 neutral, 1 cation
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          for k in 1:nLifun loop
            if datafun[k].z > 0 then
              ln_gammatau[i] :=ln_gammatau[i] + 6*m_i[j]*m_i[k]*zeta_ijktau[i, j, k];
            end if;
          end for;
        end if;
      end for;

    //neutrals
    else
      //neutral-neutral
      for j in 1:nLifun loop
        if datafun[j].z == 0 then
          ln_gammatau[i] := ln_gammatau[i] + 2*m_i[j]*lambda_ijtau[i,j];
        end if;
      end for;
      //neutral-cation
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          ln_gammatau[i] :=ln_gammatau[i] + 2*m_i[j]*lambda_ijtau[i, j];
        end if;
      end for;
      //neutral-anion
      for j in 1:nLifun loop
        if datafun[j].z < 0 then
          ln_gammatau[i] :=ln_gammatau[i] + 2*m_i[j]*lambda_ijtau[i, j];
        end if;
      end for;
      //neutral-cation-anion
      for j in 1:nLifun loop
        if datafun[j].z > 0 then
          for k in 1:nLifun loop
            if datafun[k].z < 0 then
              ln_gammatau[i] :=ln_gammatau[i] + m_i[j]*m_i[k]*zeta_ijktau[i, j, k];
            end if;
          end for;
        end if;
      end for;

    end if;
  end for;

  log_gammatau :=ln_gammatau/log(10);
end calc_loggammatau;
