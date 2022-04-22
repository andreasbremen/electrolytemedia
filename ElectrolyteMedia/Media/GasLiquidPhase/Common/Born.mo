within ElectrolyteMedia.Media.GasLiquidPhase.Common;
record Born
  constant Real a_i[:] = {14.70333593, 212.8462733, -115.4445173, 19.55210915, -83.3034798, 32.13240048, -6.694098645, -37.86202045, 68.87359646, -27.29401652};
  constant SI.Temperature T_eps = 298.15;
  constant SI.Density rho_eps = 1000;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end Born;
