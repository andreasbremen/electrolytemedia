within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solvent.tau;
function calc_log_a_ijtau "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] log_a_ijtau;

protected
  Real[nLfun-1] mol =  calc_mfromX(X);
  Real A_logtau = DebyeHueckel.Reduced.Additional.dg_dtau.calc_dAlog_dtau(T,p);
  Real I = calc_I(X);
  Real sigma=Bromley.Solvent.calc_sigma(X);

algorithm

  for i in 1:nLfun-1 loop
    for j in 1:nLfun-1 loop
      if datafun[i].z > 0 and datafun[j].z < 0 then
         log_a_ijtau[i,j] :=-1/log(10)*IF97.MH2O*(mol[i]+mol[j])*(-log(10)*(A_logtau*abs(datafun[i].z*datafun[j].z)*I^(0.5)/3*sigma));
      end if;
    end for;
  end for;

end calc_log_a_ijtau;
