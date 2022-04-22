within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Debye;
function calc_loggamma_i
  "Calculates decadic logarithm of activity coefficient in molality base of aqueous species from Debye limiting law"

  input SI.Temperature T;
  input SI.Pressure p;
  input SI.MassFraction[nLfun] X;
  output Real[nLfun-1] loggamma;
protected
  Real I=calc_I(X);
  Real A=DebyeHueckel.calc_A_log(T, p);
algorithm
  for i in 1: nLfun-1 loop
    if datafun[i].z == 0 then
      loggamma[i] := 0;//0.074 * I;
    else
      loggamma[i] := -A * datafun[i].z ^ 2 * sqrt(I);// / (1 + B * 1e-10 * 3.72 * sqrt(I));//datafun[i].a_0
    end if;
  end for;
  annotation(smoothOrder(normallyConstant = datafun)=20);
end calc_loggamma_i;
