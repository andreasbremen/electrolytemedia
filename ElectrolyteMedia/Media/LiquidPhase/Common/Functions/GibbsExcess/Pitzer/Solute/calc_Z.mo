within ElectrolyteMedia.Media.LiquidPhase.Common.Functions.GibbsExcess.Pitzer.Solute;
function calc_Z "Absolute charge molality summation"
  input Modelica.SIunits.MassFraction[nLfun] X;
  output Real Z;

protected
  Real[nLifun] m_i = calc_mfromX(X);
algorithm
  for i in 1:nLifun loop
    Z := Z + m_i [i] * abs(datafun[i].z);
  end for;
end calc_Z;

function calc_z "Unsymmetrical mixing, Pitzer p.124"
  input Real x;
  output Real z;
algorithm
  if x < 1 then
    z :=4*x^(1/5) - 2;
  else
    z :=40/9*x^(-1/10) - 22/9;
  end if;
end calc_z;
