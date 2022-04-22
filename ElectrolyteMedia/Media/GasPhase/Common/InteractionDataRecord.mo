within ElectrolyteMedia.Media.GasPhase.Common;
record InteractionDataRecord "Data record containing coefficients for interactions between Peng-Robinson gases"
  Real[:,:] k1_ij = zeros(1,1);
  Real[:,:] k2_ij = zeros(1,1);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)));
end InteractionDataRecord;
