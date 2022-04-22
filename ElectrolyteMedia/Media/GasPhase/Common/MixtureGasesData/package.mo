within ElectrolyteMedia.Media.GasPhase.Common;
package MixtureGasesData "Properties of fluids from SUPCRTBL by Zimmer(2016) and Johnson(1992)"
  constant InteractionDataRecord SingleGas(
    k1_ij={{0}},
    k2_ij={{0}}) "Single species gas";
  constant InteractionDataRecord CarbonationGas(
    k1_ij={{0,0.1896},{0.1896,0}},
    k2_ij={{0,0},{0,0}}) "Carbon dioxide,water";
  constant InteractionDataRecord H2O_CO2_O2(k1_ij={{0,0.1896,0},{0.1896,0,0},{0,0,0}},
      k2_ij={{0,0,0},{0,0,0},{0,0,0}}) "Water,carbon dioxide,oxygen";
end MixtureGasesData;
