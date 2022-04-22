within ElectrolyteMedia.Media.GasPhase.Common.Functions.PengRobinson.Molar;
function calc_a
  "Calculates product of van der Waals attraction and acentric factor in Peng Robinson EOS"

  input SI.Temperature T;
  input SI.MoleFraction[nGfun] y_i;

  output Real a( unit = "N.m4/mol2");
protected
  SI.MoleFraction[nGfun] y=y_i;
  Real[nGfun,nGfun] a_ij(unit="N.m4/mol2") = Molar.calc_a_ij(T);
  Real[nGfun,nGfun] a_ij_(unit="N.m4/mol2");
algorithm
  for i in 1:nGfun loop
    for j in 1:nGfun loop
      a_ij_[i,j] :=y[i]*y[j]*a_ij[i,j];
    end for;
  end for;
  a :=sum(a_ij_);

end calc_a;

function calc_A
  "Calculates product of van der Waals attraction and acentric factor over squared RT in Peng Robinson EOS"

  input SI.Temperature T;
  input SI.Density d;
  input SI.MoleFraction[nGfun] y_i;

  output Real A(unit="");
protected
  SI.Pressure p=Molar.calc_p(T,d,y_i);
  Real a_alpha(unit="N.m4/mol2") = Molar.calc_a(T, y_i);
algorithm
  A := a_alpha*p/((Modelica.Constants.R*T)^2);
end calc_A;
