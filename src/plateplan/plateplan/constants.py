from dataclasses import dataclass

from . import units

strains = {
    "SMU": "S. mutans UA159",
    "SGO": "S. gordonii CH1",
    "SSA": "S. sanguinus SK36",
    "SSO": "S. sobrinus SL1",
}

environments = {"AN": "anaerobic", "AE": "aerobic"}


@dataclass
class Labware:
    """
    Labware for liquid handling.
    
    id_: str
        Identifier for the labware.
    vendor: str
        Vendor name.
    catalog: str
        Catalog reference number.
    well_ids: [str]
    max_volume: units.Unit
    working_volume: units.Unit
    shape: (int, int)
    mantis_plate_filename: str
    
    """

    id_: str
    vendor: str
    catalog: str
    well_ids: [str]
    max_volume: units.Unit
    working_volume: units.Unit
    shape: (int, int)
    mantis_plate_filename: str

    def n_wells(self):
        return self.shape[0] * self.shape[1]


labware = dict()

labware["dWP96_1ml"] = Labware(
    id_="dWP96_1ml",
    vendor="VWR",
    catalog="10755-258",
    well_ids=[str(i + 1) for i in range(96)],
    max_volume=units.parse("1.2 ml"),
    working_volume=units.parse("1 ml"),
    shape=(8, 12),
    mantis_plate_filename="PT3-96-Assay-DW.pd.txt",
)

labware["dWP96_2ml"] = Labware(
    id_="dWP96_1ml",
    vendor="VWR",
    catalog="10755-248",
    well_ids=[str(i + 1) for i in range(96)],
    max_volume=units.parse("2.2 ml"),
    working_volume=units.parse("2 ml"),
    shape=(8, 12),
    mantis_plate_filename="PT3-96-Assay-DW.pd.txt",
)

labware["WP384"] = Labware(
    id_="WP384",
    vendor="Corning",
    catalog="07-201-157",
    well_ids=[str(i + 1) for i in range(384)],
    max_volume=units.parse("112 ul"),
    working_volume=units.parse("80 ul"),
    shape=(16, 24),
    mantis_plate_filename="PT9-384-Assay.pd.txt",
)

conical15 = "Conical15ml"
conical50 = "Conical50ml"
bottle250ml = "Bottle250ml"
bottle1000ml = "Bottle1000ml"

labware[conical15] = Labware(
    id_=conical15,
    vendor="Fisher",
    catalog="TBA",
    well_ids=[],
    max_volume=units.parse("16 ml"),
    working_volume=units.parse("15 ml"),
    shape=(0, 0),
    mantis_plate_filename="",
)

labware[conical50] = Labware(
    id_=conical50,
    vendor="Fisher",
    catalog="TBA",
    well_ids=[],
    max_volume=units.parse("60 ml"),
    working_volume=units.parse("50 ml"),
    shape=(0, 0),
    mantis_plate_filename="",
)

labware["Bottle250ml"] = Labware(
    id_="Bottle250ml",
    vendor="Fisher",
    catalog="TBA",
    well_ids=[],
    max_volume=units.parse("250 ml"),
    working_volume=units.parse("250 ml"),
    shape=(0, 0),
    mantis_plate_filename="",
)

labware["Bottle1000ml"] = Labware(
    id_="Bottle1000ml",
    vendor="Fisher",
    catalog="TBA",
    well_ids=[],
    max_volume=units.parse("1000 ml"),
    working_volume=units.parse("1000 ml"),
    shape=(0, 0),
    mantis_plate_filename="",
)


# LH robot configs
# {'location_id' : [allowed labware IDs]}
mantis_lc3_config = {
    "allowed": {
        "LCC(L)-1": [conical50, bottle250ml, bottle1000ml],
        "LCC(L)-2": [conical50, bottle250ml, bottle1000ml],
        "LCC(L)-3": [conical50, bottle250ml, bottle1000ml],
        "LCC(L)-4": [conical50, bottle250ml, bottle1000ml],
        "LC3(R)-1": [conical15],
        "LC3(R)-2": [conical15, conical50],
        "LC3(R)-3": [conical15],
        "LC3(R)-4": [conical15],
        "LC3(R)-5": [conical15, conical50],
        "LC3(R)-6": [conical15],
        "LC3(R)-7": [conical15],
        "LC3(R)-8": [conical15, conical50],
        "LC3(R)-9": [conical15],
        "LC3(R)-10": [conical15],
        "LC3(R)-11": [conical15, conical50],
        "LC3(R)-12": [conical15],
        "LC3(R)-13": [conical15],
        "LC3(R)-14": [conical15, conical50],
        "LC3(R)-15": [conical15],
        "LC3(R)-16": [conical15],
        "LC3(R)-17": [conical15, conical50],
        "LC3(R)-18": [conical15],
        "LC3(R)-19": [conical15],
        "LC3(R)-20": [conical15, conical50],
        "LC3(R)-21": [conical15],
        "LC3(R)-22": [conical15],
        "LC3(R)-23": [conical15, conical50],
        "LC3(R)-24": [conical15],
    },
    "priority": {
        "Conical15ml": [
            [
                "LC3(R)-1",
                "LC3(R)-4",
                "LC3(R)-7",
                "LC3(R)-10",  # Center 8 prioritized since they can
                "LC3(R)-13",
                "LC3(R)-16",
                "LC3(R)-19",
                "LC3(R)-22",
            ],  # be accessed without moving 50ml adapters
            [
                "LC3(R)-3",
                "LC3(R)-6",
                "LC3(R)-9",
                "LC3(R)-12",
                "LC3(R)-15",
                "LC3(R)-18",
                "LC3(R)-21",
                "LC3(R)-24",
            ],
        ],
        "Conical50ml": [["LCC(L)-1", "LCC(L)-2", "LCC(L)-3", "LCC(L)-4"]],
        "Bottle250ml": ["LCC(L)-4"],
        "Bottle1000ml": ["LCC(L)-4"],
    },
}
