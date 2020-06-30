"""Generate csv file mapping from well to containing solution."""

import os
import string
import csv

import pandas as pd
import numpy as np

from . import units
from . import constants


def _number_to_well(dimensions=(8, 12), reverse=False):
    num_to_well = dict()
    well_to_num = dict()

    n = 1
    for r in range(dimensions[0]):
        for c in range(1, dimensions[1] + 1):
            well = string.ascii_uppercase[r] + str(c)
            if reverse:
                well_to_num[well] = n
            else:
                num_to_well[n] = well
            n += 1

    if reverse:
        return well_to_num
    else:
        return num_to_well


def save_well_map(path, layout, instruction_set_id, dimensions):

    """Save CSV file which maps all wells to experiment performed
    including all experimental conditions.

    Inputs
    ------
    path: str, default="../data"
        Path of the parent folder of the final output file.
    layout: pd.DataFrame
        Plate layout from scheduling.schedule_mantis().
    instruction_set_id: str
        Instruction set ID from scheduling.InstructionSet().
    dimensions: tuple(int, int)
        The shape of the plate (row, col) of wells. You can use constants.shape

    Outputs
    -------
    A CSV at the location "[path]/plate_maps/map.csv" with 
    columns:
        
        parent plate: str
            Plate IDs from scheduling.LoadPlate().id.
        final plate: str
            Final plate names.
        parent_well_index: int
            1-indexed array of well IDs corresponding to the 
            parent plate.
        parent_well: str
            Well plate friendly name of parent_well_index, e.g. 5 -> A5.
        final_wells: str
            Space-separated names of well plate friendly names for the 
            corresponding wells in the final plate.   
        solution_id: str
            IDs from solution.id.
        replicate: int
            Replicate number (1-indexed).
        plate_control: bool
            True if the row is a control assay.
        plate_blank: bool
            True if the row is a blank assay.
    
    """

    num_to_well = _number_to_well(dimensions=dimensions)

    layout = layout.rename(
        columns={"well": "parent_well_index", "plate": "parent_plate"}
    )
    layout.insert(2, "parent_well", layout["parent_well_index"])
    layout["parent_well"] = layout["parent_well"].apply(lambda x: num_to_well[x])

    layout = layout.sort_values(by=["parent_plate", "parent_well_index"])

    output_folder = os.path.join(path, "plate_maps", instruction_set_id)
    output_path = os.path.join(output_folder, "map.csv")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    layout.to_csv(output_path, index=False)


def files_to_layout(filepath):
    """Converts from Mantis worklists to proper well map.
    
    WARNING: Not supported. This can be used to recover an equivalent mapping if you forget
    to use mantis.save_well_map() when running your protocol.
    
    Customize values below before using.

    """

    ####### CUSTOMIZE FILENAMES ACCORDINGLY:
    plate_dispense_files = {
        "Plate 1$fa8941a6392b8cbce99989927365caae983ba598": ["ad39d51d_mantis"],
        "Plate 2$778cc3206be303aa14a495f9e2d4dba38cf1b97d": ["cc6af787_mantis"],
        "Plate 3$0c5bfb6c6a1c28336201f450ff16b8393ad29399": ["489a3d0d_mantis"],
        "Plate 4$95a3d25c285c527f253cc3afaf1e8806e834a160": ["a32ff7ac_mantis"],
    }
    #######

    all_rows = list()
    for plate, files in plate_dispense_files.items():
        L2O_reagents = dict()
        for f in files:
            filename = os.path.join(filepath, f + ".dl.txt")
            with open(filename, newline="") as csvfile:
                reader = csv.reader(csvfile, delimiter="\t")
                no_header = list(reader)[5:]
                reagents = [no_header[idx][0] for idx in range(0, len(no_header), 10)]
                del_indexes = list()
                for idx in range(0, len(no_header), 10):
                    del_indexes.append(idx)
                    del_indexes.append(idx + 1)
                for idx in sorted(del_indexes, reverse=True):
                    del no_header[idx]

                f = pd.DataFrame(no_header).to_numpy(dtype=np.float32)
                f = np.reshape(f, (-1, 8, 12))
                result = np.where(f == 0)

                for coord in list(zip(result[0], result[1], result[2])):
                    reagent = reagents[coord[0]]
                    reagent_short = reagent[:-47] if reagent != "water" else reagent
                    well_num = coord[1] * 12 + coord[2] + 1
                    if L2O_reagents.get(well_num, False):
                        L2O_reagents[well_num].add(reagent_short)
                    else:
                        L2O_reagents[well_num] = set([reagent_short])
        for well, reagents in L2O_reagents.items():
            solution = "CDM --" + " --".join(reagents)
            new_row = [
                plate,
                well,
                solution,
                False,
                False,
                False,
                "S. mutans UA159", #### CUSTOMIZE
                "5% CO2",  #### CUSTOMIZE
            ]
            all_rows.append(new_row)

    layout = pd.DataFrame(
        all_rows,
        columns=[
            "plate",
            "well",
            "solution_id",
            "replicate",
            "plate_control",
            "plate_blank",
            "strain",
            "environment",
        ],
    )
    #### YOU'LL NEED TO SAVE THIS TO A MAP.CSV
    return layout


def collect_data(map_csv, plate_indexes, initial_data, final_data, output_path):
    map_csv = pd.read_csv(map_csv).sort_values(by=["solution_id", "replicate"])
    map_csv.insert(0, "initial_od", map_csv["parent_well_index"])
    map_csv.insert(1, "final_od", map_csv["parent_well_index"])

    initial = pd.read_excel(initial_data, sheet_name=None)
    initial_dataframes = list()
    for key, sheet in initial.items():
        initial_dataframes.append(
            sheet.drop(labels=range(0, 23))
            .drop(labels=["Unnamed: 0", "Unnamed: 1", "Unnamed: 26"], axis=1)
            .to_numpy()
            .flatten()
        )
    final = pd.read_excel(final_data, sheet_name=None)
    final_dataframes = list()
    for key, sheet in final.items():
        final_dataframes.append(
            sheet.drop(labels=range(0, 23))
            .drop(labels=["Unnamed: 0", "Unnamed: 1", "Unnamed: 26"], axis=1)
            .to_numpy()
            .flatten()
        )

    def _get_initial_od(parent_plate, parent_well_index):
        return initial_dataframes[plate_indexes[parent_plate]][parent_well_index - 1]

    def _get_final_od(parent_plate, parent_well_index):
        return final_dataframes[plate_indexes[parent_plate]][parent_well_index - 1]

    map_csv["initial_od"] = map_csv.apply(
        lambda x: _get_initial_od(x.parent_plate, x.parent_well_index), axis=1
    )
    map_csv["final_od"] = map_csv.apply(
        lambda x: _get_final_od(x.parent_plate, x.parent_well_index), axis=1
    )
    map_csv["solution_id"] = map_csv["solution_id"].apply(lambda x: x.split(" -- "))

    max_num = max(map_csv["solution_id"].apply(lambda x: len(x)))
    map_csv.insert(0, "base_solution", map_csv["solution_id"])
    map_csv["base_solution"] = map_csv["base_solution"].apply(lambda x: x[0])
    for i in range(1, max_num):
        map_csv.insert(i, f"leave_out_{i}", map_csv["solution_id"])
        map_csv[f"leave_out_{i}"] = map_csv[f"leave_out_{i}"].apply(
            lambda x: x[i] if i < len(x) else None
        )

    map_csv["final_od"] = map_csv.apply(
        lambda x: _get_final_od(x.parent_plate, x.parent_well_index), axis=1
    )

    map_csv = map_csv.drop(labels=["solution_id"], axis=1)
    map_csv.to_csv(output_path, index=False)
    return True


if __name__ == "__main__":

    save_well_map(
        layout,
        "SAVED_MAP.csv",
        num_split_plates=0,
        ot_split_file="96_to_384_Mapping_8plate.csv",
    )