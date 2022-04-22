within ElectrolyteMedia.Media.LiquidPhase.Common;
record StdRefH2O "Record containing data from Helgeson 1974"
  constant SI.SpecificEnthalpy h_tr = -287720.759/0.0180158 "data from Helgeson1974";
  constant SI.SpecificEnergy g_tr = -235512.710/0.0180158 "data from Helgeson1974";
  constant SI.SpecificEnergy u_tr = -284027.643/0.0180158 "data from Helgeson1974";
  constant SI.SpecificEnergy a_tr = -231855.624/0.0180158 "data from Helgeson1974";
  constant SI.SpecificEntropy s_tr = 63.31262/0.0180158 "data from Helgeson1974";
  constant SI.Temperature T_tr = 273.16 "triple point temperature";
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end StdRefH2O;
