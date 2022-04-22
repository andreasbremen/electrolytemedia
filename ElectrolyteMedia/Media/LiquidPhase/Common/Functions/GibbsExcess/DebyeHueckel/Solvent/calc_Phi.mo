within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.DebyeHueckel.Solvent;
function calc_Phi
  "calculates osmotic coefficient based on extended Debye-Hückel model with integration taken from Robinson and Stokes (1959)"
  input Modelica.SIunits.Temperature T;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real Phi;
protected
  Real[nLfun-1] mol = calc_mfromX(X);
  Real I=calc_I(X);
  Real A=calc_A_log(T, p);
  Real B=calc_B_log(T, p);
  Real a_0;
  Real sigma;
algorithm
  if I > 0 then
    a_0 :=sum(0.5*(datafun[i].z)^2*mol[i]*datafun[i].a_0 for i in 1:nLfun-1)/I;
    sigma :=calc_sigma(B*1e-10*a_0*sqrt(I));
  else
    sigma :=1;
  end if;

  Phi :=1 - A*sqrt(I)/3*sigma;
end calc_Phi;
