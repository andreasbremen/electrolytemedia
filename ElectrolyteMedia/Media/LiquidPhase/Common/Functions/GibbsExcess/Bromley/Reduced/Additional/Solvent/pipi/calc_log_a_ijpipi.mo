within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Bromley.Reduced.Additional.Solvent.pipi;
function calc_log_a_ijpipi "Helper function to calculate Gibbs derivative"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;

  output Real[nLfun-1,nLfun-1] log_a_ijpipi;

protected
  Real[nLfun-1] mol = calc_mfromX(X);
  Real A_logpipi = DebyeHueckel.Reduced.Additional.d2g_dpi2.calc_d2Alog_dpi2(T,p);
  Real I = calc_I(X);
  Real sigma=Bromley.Solvent.calc_sigma(X);

algorithm

  for i in 1:nLfun-1 loop
    for j in 1:nLfun-1 loop
      if datafun[i].z > 0 and datafun[j].z < 0 then
                log_a_ijpipi[i,j] :=-1/log(10)*IF97.MH2O*(mol[i]+mol[j])*(-log(10)*(A_logpipi*abs(datafun[i].z*datafun[j].z)*I^(0.5)/3*sigma));
      end if;
    end for;
  end for;

end calc_log_a_ijpipi;
