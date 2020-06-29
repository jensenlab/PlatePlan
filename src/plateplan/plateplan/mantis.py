"""Interface/file generator for the Formulatrix Mantis liquid handler

Follows the Reagent Name Based format:
    http://formulatrix.com/online_helps/mantis/4.0/Content/
        Designing_Dispenses/Tutorials/Importing_Dispense_Data.htm

"""

import shutil
import os
import csv

import pandas as pd
import numpy as np

from . import scheduling
from . import units
from . import constants


def get_dispense_order(worklist, plate_type):
    """Get dispense order of worklist from pipetting instructions

    Return a list of stock_ids in order of dispense

    worklist: pd.DataFrame
        Single plate pipetting instructions
    plate_type: constants.Labware()
        Type of plate
    """
    pivot_worklist = get_dispense_list(worklist, plate_type)
    return pivot_worklist.columns


def get_dispense_volumes(worklists, plate_type):
    """Get total volume of each stock from multiple pipetting 
       instructions

    Return a dictionary from stock_id->volume

    worklist_ids: [str]
        List of worklist CSVs to compute total volume from
    instruction_set_id: str
        ID for instruction set
    plate_type: constants.Labware()
        Type of plate
    """
    worklists = [get_dispense_list(w, plate_type) for w in worklists]
    combined = pd.concat(worklists)
    volumes = combined.sum(axis=0).round(decimals=1)
    return dict(zip(volumes.index.values, volumes.values))


def get_dispense_list(worklist, plate_type):
    """Get worklist

    Return a pd.DataFrame -> wells as index, columns as stocks, 
        and values as volume.

    worklist: pd.DataFrame
        Single plate pipetting instructions
    plate_type: constants.Labware()
        Type of plate
    """

    # Pivot with wells as index, columns as stocks, and values as volume
    pivot_worklist = worklist.pivot(index="well", columns="stock", values="volume")

    # Add rows for wells that don't appear
    indexes_present = set(pivot_worklist.index.values)
    indexes_all = set(range(1, plate_type.n_wells() + 1))
    indexes_not_present = list(indexes_all - indexes_present)
    for i in indexes_not_present:
        pivot_worklist = pivot_worklist.append(pd.Series(name=i))

    # Fill empty cells with 0 volume
    pivot_worklist = pivot_worklist.fillna(0)

    # Convert all Unit objects to strings
    def _get_value(cell):
        if isinstance(cell, units.Unit):
            return cell.convert("ul").value
        return cell

    pivot_worklist = pivot_worklist.applymap(_get_value)

    # Ensure sorted by ascending well; round decimals
    pivot_worklist = pivot_worklist.sort_index().round(decimals=1)
    return pivot_worklist


def format_worklist_dl_txt(worklist, plate_type, path):
    """Saves worklist into a Mantis-compatible protocol. Forwards compatible with
    Formulatrix Mantis software supporting Mantis worklist Version 5 or greater. Not
    backwards compatible for any software that only supports up to Version 4.

    worklist: pd.DataFrame
        Single plate pipetting instructions
    plate_type: constants.Labware()
        Type of plate
    path: str
        Output path including [filename].dl.txt
    """
    worklist = get_dispense_list(worklist, plate_type)
    reagents = worklist.columns
    num_reagents = len(reagents)
    delay_header = [num_reagents]
    for i in range(num_reagents):
        delay_header.extend([0, ""])

    with open(path, "w", newline="", encoding="utf-8") as f:
        # Header
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["[ Version: 5 ]"])
        writer.writerow([plate_type.mantis_plate_filename])
        writer.writerow(delay_header)
        writer.writerow([1])
        writer.writerow(delay_header)

        # Pipette volumes
        for idx, r in enumerate(reagents):
            writer.writerow([r, "", "Normal"])
            writer.writerow(["Well", 1])
            vols = worklist.iloc[:, idx].values
            vols = vols.reshape(plate_type.shape)
            writer.writerows(vols)


def worklist_to_csv(worklist, path):
    """Convert worklist and output a CSV at path.

    Return a list of reagents in the order they were 
    mapped from (reagent_1, reagent_N)->(1, N).

    worklist: pd.DataFrame
        Wells as index, columns as stocks, and values as volume.
    path: str
        Output path including [filename].csv
    """
    reagent_order = worklist.columns
    num_columns = len(reagent_order)

    # Removes index column, converts reagent names -> range of numbers
    # header = []
    # for i in range(0,num_columns):
    #     header.append("{} - {}".format(worklist.columns[i], i+1))
    header = pd.MultiIndex.from_arrays(
        [worklist.columns.tolist(), list(range(1, num_columns + 1))]
    )
    worklist.columns = header
    worklist.to_csv(path, index=False)
    return reagent_order


def _make_archive(source, destination):
    """
    Makes a .zip file from a source folder, and saves it in destination.
    """
    base = os.path.basename(destination)
    name = base.split(".")[0]
    format = base.split(".")[1]
    archive_from = os.path.dirname(source)
    archive_to = os.path.basename(source.strip(os.sep))
    shutil.make_archive(name, format, archive_from, archive_to)
    shutil.move("{}.{}".format(name, format), destination)


def generate_experiment_files(
    instruction_set, layout, inoculation_volume, path_to_worklists=""
):
    """Create worklists for reagent and inculation dispenses
       and save CSV files

    Save all worklists from experiment as .csv files in 
        data/worklists/[experiment ID].
    Saves a single .zip in 
        data/worklists/[experiment ID]/[experiment ID]_all.zip.

    instruction_set: scheduling.InstructionSet
        Instruction set for the experiment.
    layout: pd.Dataframe
        Plate layout, mapping every well to its content for the 
        experiment.
    innoculation_volume: str
        To be interpreted with unit.parse() to determine inoculation
        volume.
    path_to_worklists: str, default = ""
        Path to locate the 'worklist' folder where files will be saved.
        Default is current directory.
    """
    _generate_reagent_worklists(instruction_set, path_to_worklists)
    _generate_inoculation_worklists(
        instruction_set, layout, inoculation_volume, path_to_worklists
    )

    path = os.path.join(path_to_worklists, "worklists")
    new_dir = os.path.join(path, instruction_set.id)
    archive_dest = os.path.join(new_dir, "{}_all.zip".format(instruction_set.id))
    if not os.path.exists(path):
        os.mkdir(path)
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    _make_archive(new_dir, archive_dest)


def _generate_reagent_worklists(instruction_set, path_to_worklists=""):
    """Convert worklist and save CSV files

    Save all worklists from experiment as .csv files in 
        data/worklists/[experiment ID].
    Saves a single .zip in 
        data/worklists/[experiment ID]/[experiment ID]_all.zip.

    instruction_set: scheduling.InstructionSet
        Instruction set for the experiment.
    path_to_worklists: str, default = ""
        Path to locate the 'worklist' folder where files will be saved.
        Default is current directory.
    """

    print(("Exporting .csv and .zip files for {}".format(instruction_set.id)))

    path = os.path.join(path_to_worklists, "worklists")
    instructions = instruction_set.instructions
    new_dir = os.path.join(path, instruction_set.id)
    if not os.path.exists(path):
        os.mkdir(path)
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    for step in instructions:
        if isinstance(step, scheduling.RunMantis):
            for plate in step.plates:
                plate_type = constants.labware[plate[0].plate_type]
                if len(plate) > 1:
                    if isinstance(plate[1], scheduling.Pipette):
                        file_id = plate[1].file_id
                        worklist = plate[1].worklist
                        text_dest = os.path.join(
                            new_dir, "{}_mantis.dl.txt".format(file_id)
                        )
                        format_worklist_dl_txt(worklist, plate_type, text_dest)
                        worklist = get_dispense_list(worklist, plate_type)
                        if not os.path.exists(new_dir):
                            os.mkdir(new_dir)
                        csv_dest = os.path.join(new_dir, "{}.csv".format(file_id))
                        reagents_in_order = worklist_to_csv(worklist, csv_dest)


def _generate_inoculation_worklists(
    instruction_set, layout, inoculation_volume, path_to_worklists=""
):
    """Generate inculation dispense files from layout and save CSV files
    
    instruction_set: scheduling.InstructionSet
        Instruction set for experiment.
    layout: pd.Dataframe
        Plate layout, mapping every well to its content for the 
        experiment.
    innoculation_volume: str
        To be interpreted with unit.parse() to determine inoculation
        volume.
    path_to_worklists: str, default = ""
        Path to locate the 'worklist' folder where files will be saved.
        Default is current directory.
    """

    path = os.path.join(path_to_worklists, "worklists")
    instructions = instruction_set.instructions
    new_dir = os.path.join(path, instruction_set.id)
    if not os.path.exists(path):
        os.mkdir(path)
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)

    plate_names = {p: None for p in set(layout["plate"].values)}
    for instruction in instruction_set.instructions:
        if isinstance(instruction, scheduling.RunMantis):
            for plate in instruction.plates:
                plate_names[plate[0].plate_id] = plate[0].plate_type

    for name, plate_type in plate_names.items():
        plate_type = constants.labware[plate_type]
        plate_layout = layout[layout["plate"] == name]

        actual_wells = set(plate_layout["well"].tolist())
        all_wells = set(range(1, plate_type.n_wells() + 1))
        empty_wells = list(all_wells - actual_wells)
        worklist = plate_layout[plate_layout["plate_blank"] == True]
        zero_wells = worklist["well"].tolist() + empty_wells

        new_rows = pd.DataFrame(
            np.ones((plate_type.n_wells(), 3)), columns=["well", "stock", "volume"]
        )
        new_rows["well"] = all_wells
        new_rows["stock"] = worklist.iloc[0, 6]
        new_rows["volume"] = units.parse(inoculation_volume)
        new_rows.loc[new_rows["well"].isin(zero_wells), ["volume"]] = units.parse(
            "0 ul"
        )

        text_dest = os.path.join(new_dir, "{}_mantis.dl.txt".format(name))
        format_worklist_dl_txt(new_rows, plate_type, text_dest)


if __name__ == "__main__":

    # worklist = pd.read_pickle('worklist.pkl')
    # worklist = get_dispense_list(worklist, constants.labware['WP384'])
    # reagents_in_order = worklist_to_csv(worklist, 'worklist.csv')
    # print(worklist)
    # print(reagents_in_order)
    import protocols.CDM as CDM
    import makeids

    plates, instructions, layout = CDM.schedule_CDM_l2o(
        CDM.CDM, list(CDM.CDM_stocks.values())
    )
    instruction_set = scheduling.InstructionSet(
        id_=makeids.unique_id(prefix="test"), instructions=instructions, owner="Adam"
    )
    generate_experiment_files(instruction_set, layout, "2 ul", path_to_worklists="")
