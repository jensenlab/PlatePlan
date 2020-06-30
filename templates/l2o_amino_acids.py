import csv
import logging
import datetime
import argparse
import os
import sys
import copy

import pandas as pd

from plateplan import (
    constants,
    ingredients,
    makeids,
    mantis,
    mapper,
    scheduling,
    units,
    utils,
)

parser = argparse.ArgumentParser(
    description="Generate an experiment using the JenDrop software package. (c) 2020 Jensen Lab."
)
# parser.add_argument(
#     "-db",
#     "--use_db",
#     action = 'store',
#     type = str
#     help="To use a database to read and store Ingredients, pass in a pickled databased genereted from database.Database",
# )

args = parser.parse_args()

CDM_groups = {
    "amino_acids+NH4": set(
        [
            "dl_alanine",
            "l_arginine",
            "l_aspartic_acid",
            "l_asparagine",
            "l_cystine",
            "l_cysteine",
            "l_glutamic_acid",
            "l_glutamine",
            "glycine",
            "l_histidine",
            "l_isoleucine",
            "l_leucine",
            "l_lysine",
            "l_methionine",
            "l_phenylalanine",
            "prolines",
            "l_serine",
            "l_threonine",
            "l_tryptophan",
            "l_tyrosine",
            "l_valine",
            "ammoniums",
        ]
    ),
    "vitamins": set(
        [
            "paba",
            "biotin",
            "folic_acid",
            "niacinamide",
            "nadp",
            "pantothenate",
            "pyridoxes",
            "riboflavin",
            "thiamine",
            "vitamin_b12",
            "adenine",
            "guanine",
            "uracil",
        ]
    ),
    "salts": set(
        [
            "iron_nitrate",
            "magnesium_sulfate",
            "manganese_sulfate",
            "sodium_acetate",
            "calcium_chloride",
            "sodium_bicarbonate",
            "potassium_phosphates",
            "sodium_phosphates",
        ]
    ),
}


def make_CDM():
    reagents = dict()
    stocks = dict()
    concentrations = dict()
    amino_acid_final_concentrations = dict()
    reagent_concentrations_CDM_base = dict()
    molecular_weights = dict()

    with open(
        "files/CDM_reagents.csv", newline="", encoding="utf-8-sig"
    ) as csvfile:
        reader = csv.DictReader(csvfile)

        # Go through each ingredient in the CSV file
        for row in reader:
            short_id = row["id"]

            id_ = f"{short_id}_{row['stock_concentration_x']}x"
            mw = float(row["mw_g_mol"])
            molecular_weights[id_] = mw

            concentrations[id_] = units.parse(
                row["concentration_g_l"], "g/l"
            )  # Store the final concentration for full CDM Solution
            stock_conc = units.parse(
                row["stock_concentration_g_l"], "g/l"
            )  # Store the stock concentration

            # Add non-amino acids+NH4 to CDM Base
            if short_id in CDM_groups["amino_acids+NH4"]:
                amino_acid_final_concentrations[id_] = concentrations[id_].convert(
                    "mM", molar_mass=mw
                )  # Store the final concentration for use later when rescaling
                reagent_concentrations = {id_: stock_conc}

                #### Creating individual dispense ingredients for Non-CDM base items
                # Reagent
                reagents[id_] = ingredients.Reagent(
                    id_=id_, name=row["name"], molecular_weight=mw
                )

                # Solution
                solution = ingredients.Solution(
                    reagents=reagent_concentrations, id_=id_ + "_stock"
                )

                # Stock
                stock = ingredients.Stock(
                    ingredient=solution,
                    date_made=datetime.date.today(),
                    date_expires=datetime.date.today()
                    + datetime.timedelta(6 * 365 / 12),
                    quantity=units.parse(row["stock_volume_ml"], "ml"),
                    labware=row["labware"],
                )

                stocks[stock.id] = stock

            else:
                # Store the stock concentration for the CDM base Solution
                reagent_concentrations_CDM_base[id_] = stock_conc

        #### Create CDM_base ingredient (CDM without amino acids and ammoniums)
        # Reagent
        reagents["CDM_base"] = ingredients.Reagent(
            id_="CDM_base", name=row["name"], molecular_weight=0
        )

        # Solution
        solution = ingredients.Solution(
            reagents=reagent_concentrations_CDM_base, id_="CDM_base" + "_stock"
        )

        # Stock
        cdm_stock = ingredients.Stock(
            ingredient=solution,
            date_made=datetime.date.today(),
            date_expires=datetime.date.today() + datetime.timedelta(6 * 365 / 12),
            quantity=units.parse("10 ml"),
            labware=row["labware"],
        )

        stocks[cdm_stock.id] = cdm_stock

        #### Create CDM final solution -- what we are basing our dispenses on
        cdm = ingredients.Solution(concentrations, id_="CDM")

    return reagents, stocks, cdm, amino_acid_final_concentrations, molecular_weights


def scale_NH4_concentrations(to_remove, final_conc_in, molecular_weights):
    # returns dict of molar concentrations scaled to compensate for the one removed
    final_conc = pd.DataFrame.from_dict(final_conc_in, orient="index")
    unit = final_conc.iloc[0, 0].unit_str
    final_conc[0] = final_conc[0].apply(lambda x: x.convert(unit))
    final_conc.insert(1, "value", final_conc[0])
    final_conc.insert(1, "unit", final_conc[0])
    final_conc["value"] = final_conc["value"].apply(lambda x: x.value)
    final_conc["unit"] = final_conc["unit"].apply(lambda x: x.unit_str)
    final_conc = final_conc.drop(0, axis=1)
    to_remove = utils.assert_list(to_remove)
    removed_conc = sum([final_conc.loc[r].value for r in to_remove])
    final_conc.loc["ammoniums_100x", "value"] = removed_conc / 2
    final_conc = final_conc.drop(to_remove, axis=0)
    final_conc.insert(0, "final", pd.Series(dtype=str))
    final_conc["final"] = final_conc[["value", "unit"]].apply(
        lambda x: (
            units.parse(x[0], x[1]).convert("g/l", molar_mass=molecular_weights[x.name])
        ),
        axis=1,
    )
    new_concentrations = final_conc.to_dict()["final"]
    return new_concentrations


def schedule_CDM_l2o(
    CDM, stocks, amino_acid_final_concentrations, molecular_weights,
):

    # Configure each Solution to be dispensed according to the
    # experiments we want to do.
    aas = set(amino_acid_final_concentrations.keys())
    aas.remove("ammoniums_100x")
    ids = list(aas)
    n = len(ids)
    solutions = [CDM]

    # Setting up L1O experiments for amino acids.
    # Make a solution of CDM with singular AA removed.
    # Add each to a list of solutions.
    for x in ids:
        # We decided to add a molar equivalent of ammonium solution
        # to account for the loss of nitrogen when removing an amino acid
        new_conc = scale_NH4_concentrations(
            x, amino_acid_final_concentrations, molecular_weights
        )
        new = CDM.without([x])  # Generates a new Solution object without `x` reagent
        new = new.updated(new_conc)  # updates concentrations
        solutions.append(new)

    # Setting up L2O experiments for amino acids.
    # Make a solution of CDM with every pair of AA removed.
    # Add each to a list of solutions.
    for x in range(n):
        for y in range(x + 1, n):
            # We decided to add a molar equivalent of ammonium solution
            # to account for the loss of nitrogen when removing an amino acid
            new_conc = scale_NH4_concentrations(
                [ids[x], ids[y]], amino_acid_final_concentrations, molecular_weights
            )
            new = CDM.without(
                [ids[x], ids[y]]
            )  # Generates a new Solution object without `x` and `y` reagents
            new = new.updated(new_conc)  # updates concentrations
            solutions.append(new)

    print("# of solutions:", len(solutions))

    # Actual call to the scheduling program
    return scheduling.schedule_mantis(
        solutions,
        [constants.strains["SMU"]],
        ["aerobic"],
        stocks,
        excess="water",
        plate=constants.labware["WP384"],
        total_volume="80 ul",
        working_volume="78 ul",
        plate_control=("CDM", 3),
        plate_blank=("CDM", 3),
        replicates=4,
        min_drop_size="0.1 ul",
    )


if __name__ == "__main__":
    scheduling.logger.setLevel(logging.INFO)
    (
        CDM_reagents,
        CDM_stocks,
        CDM,
        amino_acid_final_concentrations,
        molecular_weights,
    ) = make_CDM()

    # Schedule the experiment
    plates, instructions, layout = schedule_CDM_l2o(
        CDM,
        list(CDM_stocks.values()),
        amino_acid_final_concentrations,
        molecular_weights,
    )

    # Make an InstructionSet based on the instructions outputted from
    # the scheduler above. Used for generating the files below.
    friendly_name = "CDM-L2O"
    exp_id = makeids.unique_id(prefix=friendly_name)
    instruction_set = scheduling.InstructionSet(
        id_=exp_id, instructions=instructions, owner="Adam",
    )

    # Export the files needed for the experiment
    output_filepath = os.path.join("experiment_outputs", exp_id)
    os.makedirs(output_filepath)

    # Export files in Mantis-friendly format
    mantis.generate_experiment_files(
        instruction_set, layout, "2 ul", path_to_worklists=output_filepath
    )
    # Export an experiment map for future use when mapping data back to each experiment
    mapper.save_well_map(
        path=output_filepath,
        layout=layout,
        instruction_set_id=instruction_set.id,
        dimensions=constants.labware["WP384"].shape,
    ) # Saves map.csv
