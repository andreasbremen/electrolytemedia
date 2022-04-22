within ElectrolyteMedia.Media.Common.Types;
type LiquidModel = enumeration(
    Ideal         "Ideal liquid phase at infinite dilution",
    Debye         "Debye limiting law",
    DebyeHueckel  "Debye-Hueckel model",
    Bromley       "Bromley model",
    Pitzer        "Pitzer model") "Lists all types of liquid models";
