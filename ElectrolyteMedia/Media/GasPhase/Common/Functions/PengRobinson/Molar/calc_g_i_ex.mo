within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_g_i_ex
  "Calculates excess Gibbs free energy of gas species"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output SI.MolarEnergy g_ex[nGfun];
protected
  SI.Pressure[nGfun] f_i;
  SI.Pressure[nGfun] p_i;
  Real[nGfun] phi;
  SI.Pressure p;
algorithm
  p := Molar.calc_p(T,d,y_i);
  phi := Molar.calc_phi(T,d,y_i);
  p_i :=y_i*p;
  f_i :=phi .* p_i;
  for i in 1:nGfun loop
    if f_i[i] > 0 then
      g_ex[i] := Modelica.Constants.R * T * log(f_i[i]/pref);
    else
      g_ex[i] :=0;
    end if;
  end for;

end calc_g_i_ex;
