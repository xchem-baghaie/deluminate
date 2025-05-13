import numpy as np

MAX_ALLOWED_MONOSYNTHON_MW_NON_BRO5 = 400
MAX_ALLOWED_MONOSYNTHON_MW_BRO5 = 700

RDKIT_PROPERTIES = [
    "amw",
    "CrippenClogP",
    "FractionCSP3",
    "NumHBD",
    "NumHBA",
    "tpsa",
    "NumRotatableBonds",
]

PROPERTY_SYNONYMS = {
    "amw": "MW",
    "CrippenClogP": "cLogP",
    "FractionCSP3": "Fsp3",
    "NumHBD": "HBD",
    "NumHBA": "HBA",
    "tpsa": "tPSA",
    "NumRotatableBonds": "RotB",
    "CA_MW": "MW",
    "CA_cLogP": "cLogP",
    "CA_cLogD_7.4": "cLogD (pH 7.4)",
    "CA_FSP3": "Fsp3",
    "CA_HBD": "HBD",
    "CA_HBA": "HBA",
    "CA_tPSA": "tPSA",
    "CA_ROT": "RotB",
    "CA_PKA": "pKa",
    "CA_PKB": "pKb",
    "CA_RING_COUNT": "Ring Count",
    "CA_AROMATIC_RING_COUNT": "Aromatic Ring Count",
    "CA_HAC": "Heavy Atom Count",
    "CA_RULEOF5": "Lipinski Compliance",
    "CA_MPO_LUNDBECK": "Lundbeck CNS MPO Score",
    "CA_MPO_PFIZER": "Pfizer CNS MPO Score",
    "Rod_All": "Rod-Likeness",
    "Disc_All": "Disc-Likeness",
    "Sphere_All": "Sphere-Likeness",
}

PMI_XY_COLUMN_NAMES = [
    "PMI_X_Avg",
    "PMI_Y_Avg",
    "PMI_X_All",
    "PMI_Y_All",
]

PMI_RDS_COLUMN_NAMES = [
    "Rod_Avg",
    "Disc_Avg",
    "Sphere_Avg",
    "Rod_All",
    "Disc_All",
    "Sphere_All",
]

VALID_PROPERTIES = PROPERTY_SYNONYMS
for val in list(set(list(PROPERTY_SYNONYMS.values()))) + list(
    set(list(PMI_XY_COLUMN_NAMES))
):
    VALID_PROPERTIES[val] = val

PROPERTIES_FOR_GRID = [
    "CA_MW",
    "CA_cLogP",
    "CA_cLogD_7.4",
    "CA_FSP3",
    "CA_HBD",
    "CA_HBA",
    "CA_tPSA",
    "CA_ROT",
]

BIN_DCT = {
    "MW": np.arange(250, 800, 50),
    "cLogP": np.arange(-5, 10, 2),
    "cLogD (pH 7.4)": np.arange(-5, 10, 2),
    "Fsp3": [round(x, 1) for x in np.arange(0, 1.1, 0.1)],
    "HBD": np.arange(0, 8),
    "HBA": np.arange(0, 16, 2),
    "tPSA": np.arange(0, 220, 20),
    "RotB": np.arange(0, 22, 2),
    "pKa": [x for x in np.arange(-6.0, 21.0, 3.0)],
    "pKb": [x for x in np.arange(-6.0, 21.0, 3.0)],
    "Ring Count": np.arange(0, 11),
    "Aromatic Ring Count": np.arange(0, 11),
    "Heavy Atom Count": np.arange(10, 90, 10),
    "Lundbeck CNS MPO Score": [round(x, 1) for x in np.arange(0, 6.5, 0.5)],
    "Pfizer CNS MPO Score": [round(x, 1) for x in np.arange(0, 6.5, 0.5)],
    "Rod-Likeness": [round(x, 1) for x in np.arange(0, 1.1, 0.1)],
    "Disc-Likeness": [round(x, 1) for x in np.arange(0, 1.1, 0.1)],
    "Sphere-Likeness": [round(x, 1) for x in np.arange(0, 1.1, 0.1)],
}

COLOR_DCT = {
    "MW": "#CC99FF",
    "cLogP": "#5C89E4",
    "cLogD (pH 7.4)": "#AFD5AF",
    "Fsp3": "#9E389E",
    "HBD": "#03A32B",
    "HBA": "#6FDC6F",
    "tPSA": "#265F97",
    "RotB": "#67B1F9",
    "pKa": None,
    "pKb": None,
    "Ring Count": None,
    "Aromatic Ring Count": None,
    "Heavy Atom Count": None,
    "Lundbeck CNS MPO Score": None,
    "Pfizer CNS MPO Score": None,
    "Rod-Likeness": None,
    "Disc-Likeness": None,
    "Sphere-Likeness": None,
}

ROUNDING_DCT = {
    "MW": 0,
    "cLogP": 2,
    "cLogD (pH 7.4)": 2,
    "Fsp3": 2,
    "HBD": 0,
    "HBA": 0,
    "tPSA": 0,
    "RotB": 0,
    "pKa": 2,
    "pKb": 2,
    "Ring Count": 0,
    "Aromatic Ring Count": 0,
    "Heavy Atom Count": 0,
}
