within ElectrolyteMedia.Media.GasLiquidPhase.Common;
record IF97 "record for IF97 parameters of region 1; pT explicit"
  extends Modelica.Media.Water.IF97_Utilities.BaseIF97.data(MH2O = 0.018016);
  constant SI.Temperature Ttriple=273.16 "The triple point temperature";
  constant SI.Pressure ptriple=611.657 "The triple point pressure";
  constant SI.Density dltriple=999.792520031617642
    "The triple point liquid density";
  constant SI.Density dvtriple=0.485457572477861372e-2
    "The triple point vapour density";
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end IF97;
