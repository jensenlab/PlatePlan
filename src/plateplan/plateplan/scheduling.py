"""Schedule robots, instruments, and assays."""

import copy
import logging
import random

import gurobipy as grb
import pandas as pd
import numpy as np

from . import makeids
from . import units
from . import utils
from . import constants


logger = logging.getLogger(__name__)
_INFO = logger.info
_DEBUG = logger.debug

USE_MULTIPROCESSING = True

INDENT = "    "


def get_indent(n):
    return "".join([INDENT for _ in range(n)])


class SchedulingError(Exception):
    def __init__(self, message):
        self.message = message


class Instruction(object):
    pass


class InstructionSet(Instruction):
    def __init__(self, id_, instructions, owner):
        self.id = id_
        self.instructions = instructions
        self.owner = owner

    def __repr__(self):
        return "InstructionSet({} stock sets)".format(len(self.instructions))

    def report(self, indent=0):
        pre = get_indent(indent)
        print(pre + "Instruction Set " + makeids.trim_id(self.id, dots=True))
        for instruction in self.instructions:
            instruction.report(indent + 1)


class LoadPlate(Instruction):
    def __init__(self, plate_id, plate_type):
        self.plate_id = plate_id
        self.plate_type = plate_type

    def __repr__(self):
        return "LoadPlate({}: {})".format(self.plate_id, self.plate_type)


# =========================== Mantis instructions ===========================


class LoadStocks(Instruction):
    def __init__(self, stock_ids):
        self.stock_ids = stock_ids

    def __repr__(self):
        return "LoadStocks({} stocks)".format(len(self.stock_ids))

    def report(self, indent=0):
        pre = get_indent(indent)
        print(pre + "Load Stocks:")
        for id_ in self.stock_ids:
            print(pre + INDENT + makeids.trim_id(id_, dots=True))


class Pipette(Instruction):
    def __init__(self, worklist, file_id=""):
        self.worklist = worklist
        self.file_id = file_id

    def __repr__(self):
        return "Pipette({} instructions)".format(len(self.worklist))

    def report(self, indent=0):
        pre = get_indent(indent)
        print(pre + "{} pipetting commands.".format(len(self.worklist)))


class RunMantis(Instruction):
    """Pipette using the Mantis liquid handler.
        
    stock: LoadStocks
        A single LoadStocks object.
    plates: [(LoadPlate, Pipette), ...]
        A list of 2-tuples with a LoadPlate object and a Pipette object.
    
    """

    def __init__(self, stock, plates):
        self.stock = stock
        self.plates = plates

    def filter_by_plates(plates):
        pass

    def report(self, indent=0):
        pre = get_indent(indent)
        print(pre + "========== Run Mantis ==========")
        self.stock.report(indent + 1)
        for plate, instructions in self.plates:
            plate_id = makeids.trim_id(plate.plate_id, dots=True)
            print(pre + INDENT + "Pipette into " + plate_id)
            instructions.report(indent + 2)


# =========================== scheduling routines ===========================


# ===========================================================================
# ========================== pipetting strategies ===========================
# ===========================================================================


def _make_mantis_pipetting_strategy(
    target_quantities,
    stock_quantities_per_drop,
    drop_size,
    working_volume,
    ingredient_ids,
    stock_ids,
    solution_id,
    require_nonzero=True,
    minimize_stocks=True,
    opt_frac=0.01,
    excess=None,
    quiet=True,
):
    """Find a Mantis pipetting strategy for a single experiment.

    Inputs
    ------
    target_quantity: np.array len=n_ingredients
        Required quantity of each reagent in the final mixture.
    quantity_per_drop: np.array dim=(n_ingredients, n_stocks)
        Quantity added per drop from each solution.
    dropsize: float
        Drop volume in ul.
    working_volume: float
        Maximum well volume.
    ingredient_ids: List[str]
        Unique IDs for each ingredient.
    stock_ids: List[str]
        Unique IDs for each stock.
    solution_id: str
        Unique ID for the solution.
    require_nonzero: bool, default=True
        If true, every ingredient must be added. This constraint applies in
        two situations:
            1. the concentration of an ingredient's stock is too high, so the
            best solution would be to add none.
            2. the concentrations of some solutions are so low that the
            maximum volume constraint is hit so some ingredients are left out.
    minimize_stocks: bool, default=True
        If true, use the fewest number of stocks while keeping the accuracy
        within `opt_frac`.
    opt_frac: float, default=0.01
        The accuracy of the minimized pipetting solution must be within
        this fraction of the solution using any number of stocks.
    excess: str, default=None
        Stock ID for the solution used to "top off" each well to reach the
        working volume. If None, the wells may have non-uniform volumes.
    quiet: boolean, default=False
        If true, no output from the MILP solver is printed.

    Returns
    -------
    A pd.DataFrame with one row and the following columns:
        quantity_data: pd.DataFrame
            Row index are ingredient IDs. Columns are
                quantity_added: total quantity added to the well
                target_quantity: target quantity specified in the inputs
                percent_error: |added - target| / target * 100
        volume_added: pd.Series[float]
            Volume of stock added to the well. Row index are stock IDs. Only
            stocks used (those with nonzero volumes added) are returned. Each
            entry is a Unit object.
    The dataframe uses the `solution_id` as the row index.
    """
    allowed_ingredients = (target_quantities > 0).nonzero()[0]
    # select only stocks with all ingredients in this assay
    stock_has_ingredient = (stock_quantities_per_drop > 0).astype(int)
    allowed_stocks = (
        stock_has_ingredient @ (target_quantities > 0)
        == stock_has_ingredient.sum(axis=1)
    ).nonzero()[0]

    target_quantities = target_quantities[allowed_ingredients]
    quantity_per_drop = (
        stock_quantities_per_drop.take(allowed_stocks, axis=0)
        .take(allowed_ingredients, axis=1)
        .transpose()
    )
    n_ingredients, n_stocks = quantity_per_drop.shape
    volume_per_drop = drop_size * np.ones(len(allowed_stocks))

    model = grb.Model()
    ub = {i: working_volume // volume_per_drop[i] for i in range(n_stocks)}
    vars = model.addVars(n_stocks, vtype=grb.GRB.INTEGER, lb=0, ub=ub)
    slacks = model.addVars(
        n_ingredients,
        vtype=grb.GRB.CONTINUOUS,
        lb=-grb.GRB.INFINITY,
        ub=grb.GRB.INFINITY,
    )
    # Sum_{stocks}(quantity_per_drop*number_of_drops) - slack = target_quantity
    for ingredient in range(n_ingredients):
        linex = grb.LinExpr()
        linex.addTerms(quantity_per_drop[ingredient, :].tolist(), vars.values())
        linex -= slacks[ingredient]
        model.addLConstr(linex, sense=grb.GRB.EQUAL, rhs=target_quantities[ingredient])
        if require_nonzero:
            # Sum_{stocks w/ ingredient}(number_of_drop) >= 1
            linex = grb.LinExpr()
            linex.addTerms(
                (quantity_per_drop[ingredient, :] > 0).astype(int).tolist(),
                vars.values(),
            )
            model.addLConstr(linex, sense=grb.GRB.GREATER_EQUAL, rhs=1)
    # Sum_{stocks}(volume_per_drop*number_of_drops) <= max_volume
    model.addLConstr(
        grb.LinExpr(volume_per_drop, vars.values()),
        sense=grb.GRB.LESS_EQUAL,
        rhs=working_volume,
    )
    # min Sum_{slacks}(weight*slack^2)
    #   = min Sum_{ingredients}(relative_error(ingredient, target)^2)
    qobj = grb.QuadExpr()
    weights = [1 / target_quantities[i] ** 2 for i in range(n_ingredients)]
    qobj.addTerms(weights, slacks.values(), slacks.values())
    model.setObjective(qobj, grb.GRB.MINIMIZE)

    if not quiet:
        model.update()
        model.printStats()
    if quiet:
        model.setParam("OutputFlag", 0)

    model.optimize()
    if model.getAttr("Status") != grb.GRB.OPTIMAL:
        if not quiet:
            model.write("out.lp")
        raise SchedulingError(
            "no pipetting strategy found; " "problem written to out.lp"
        )

    if minimize_stocks:
        # add binary indicators for each stock
        indicators = model.addVars(n_stocks, vtype=grb.GRB.BINARY)
        for stock in range(n_stocks):
            model.addGenConstrIndicator(
                indicators[stock], False, vars[stock], grb.GRB.EQUAL, 0.0
            )
        # set the slacks to keep the solution with the tolerance
        for ingredient in range(n_ingredients):
            delta = opt_frac * target_quantities[ingredient]
            current_slack = abs(slacks[ingredient].getAttr("x"))
            slacks[ingredient].setAttr("lb", -current_slack - delta)
            slacks[ingredient].setAttr("ub", current_slack + delta)
        model.setObjective(grb.quicksum(indicators), grb.GRB.MINIMIZE)
        model.optimize()
        if not quiet:
            model.write("out2.lp")
        if model.getAttr("Status") == grb.GRB.OPTIMAL:
            for indicator in indicators.values():
                indicator.setAttr("lb", np.round(indicator.getAttr("x")))
                indicator.setAttr("ub", np.round(indicator.getAttr("x")))
        model.setObjective(qobj, grb.GRB.MINIMIZE)
        if not quiet:
            model.write("out3.lp")
        model.optimize()
        if model.getAttr("Status") != grb.GRB.OPTIMAL:
            raise SchedulingError(
                "no pipetting strategy found; " "after minimizing stock count"
            )

    drops = np.around(np.array([v.getAttr("x") for v in vars.values()]))
    quantity_added = quantity_per_drop @ drops
    percent_error = np.abs(quantity_added - target_quantities) / target_quantities * 100
    quantities = pd.DataFrame(
        {
            "quantity_added": quantity_added,
            "target_quantity": target_quantities,
            "percent_error": percent_error,
        },
        index=[ingredient_ids[k] for k in allowed_ingredients],
    )

    volume_added = pd.Series(
        drops * np.array(volume_per_drop), index=[stock_ids[k] for k in allowed_stocks]
    )
    volume_added = volume_added[volume_added > 0]

    if excess is not None:
        vol = sum(volume_added)
        if working_volume - vol >= drop_size:
            to_add = drop_size * ((working_volume - vol) // drop_size)
            if excess not in volume_added:
                volume_added[excess] = to_add
            else:
                volume_added[excess] += to_add

    # convert to Units
    volume_added = pd.Series(
        [units.parse(vol, "ul") for vol in volume_added], index=volume_added.index
    )

    # I can't find a way to create a one-row DataFrame with a nested Series
    # or DataFrame. All I can do is create a 2-row DataFrame and grab the
    # first row. #ThanksPandas
    df = pd.DataFrame(
        dict(
            volume_added=[volume_added, volume_added],
            quantities=[quantities, quantities],
        ),
        index=[solution_id, "<null>"],
    )
    return df.iloc[0:1, :]


def plan_mantis_pipetting(
    solutions,
    stocks,
    total_volume="200 ul",
    working_volume="200 ul",
    excess=None,
    min_drop_size="0.1 ul",
    quantity_units="ug/ul",
):

    # TODO: check that the volumes and drop sizes have the same units
    total_volume = units.parse(total_volume).convert("ul")
    working_volume = units.parse(working_volume).convert("ul")
    min_drop_size = units.parse(min_drop_size).convert("ul")

    n_solutions = len(solutions)
    ingredient_ids = set()
    for solution in solutions:
        ingredient_ids |= solution.reagents.keys()
    n_ingredients = len(ingredient_ids)

    # remove any stock that contains an ingredient that appears in no assay
    stocks = [s for s in stocks if not set(s.ingredient_ids) - ingredient_ids]
    n_stocks = len(stocks)

    iidx = {ing: k for k, ing in enumerate(ingredient_ids)}
    # asolution_quantities[i,j] is the quantity of ingredient j in assay i
    solution_quantities = np.zeros((n_solutions, n_ingredients))
    for i, solution in enumerate(solutions):
        for ing, conc in solution.reagents.items():
            solution_quantities[i, iidx[ing]] = (
                total_volume * conc.convert(quantity_units)
            ).value

    stock_quantities_per_drop = np.zeros((n_stocks, n_ingredients))
    for i, stock in enumerate(stocks):
        for ing, conc in stock.ingredient.reagents.items():
            stock_quantities_per_drop[i, iidx[ing]] = (
                min_drop_size * conc.convert(quantity_units)
            ).value

    solution_ids = [solution.id for solution in solutions]
    stock_ids = [stock.id for stock in stocks]
    ingredient_ids = list(ingredient_ids)

    # Todo: add multiprocessing here
    dfs = []
    for k in range(n_solutions):
        dfs.append(
            _make_mantis_pipetting_strategy(
                solution_quantities[k, :],
                stock_quantities_per_drop,
                min_drop_size.value,
                working_volume.value,
                ingredient_ids,
                stock_ids,
                solution_ids[k],
                excess=excess,
            )
        )

    return pd.concat(dfs)


# ===========================================================================
# ============================== plate layouts ==============================
# ===========================================================================


def _make_locations(
    n_assays,
    wells_per_plate,
    replicates=1,
    n_plate_controls=1,
    n_plate_blanks=1,
    randomize=True,
):
    """Make lists of locations for controls, assays, and blanks.

    Inputs
    ------
    n_assays: int
        Total number of assays. This includes the control (if plate_control =
        True) but does not include replicates.
    wells_per_plate: int
        Number of wells per plate.
    replicates: int, default=1
        Number of times each assay will be assigned a location. Replicates for
        a single assay appear on the sample plate. If `randomize` is False,
        the replicates appear in adjacent wells; otherwise, the replicates
        are scattered across the plate.
    n_plate_controls: int, default=1
        If > 0, one of the assays is a control that should be repeated 
        `n_plate_controls` times on every plate.
    n_plate_blanks: int, default=1
        If > 0, one of the assays is a blank that should be repeated 
        `n_plate_blanks` times on every plate.
    randomize: bool, default=True
        If True, the assays and replicates are assigned to random wells on the
        plate, although each assay's replicates remain on the same plate.

    Returns
    -------
    A dict with the number of plates (`n_plates`) and locations for the
    assays, controls, blanks, and empty wells.
        - Plates and wells start at index 1.
        - Each location is a 2-tuple with the plate index and a list of the
        well indexes for the replicates.

    Example
    -------
    >>> _make_locations(
            wells_per_plate = 96,
            n_assays = 25,
            replicates = 4,
            plate_control = True,
            plate_blank = False)

    With replicates, this set of assays uses 25 * 4 = 100 wells. This requires
    two 96 well plates. Since a plate control is included, an additional four
    wells will be used to duplicate the control on the second plate. Thus the
    entire set of assays uses 104 wells. These are split evenly across both
    plates (52 in each).
    """

    assays_per_plate = wells_per_plate // replicates
    if n_plate_controls:
        assays_per_plate -= n_plate_controls
        n_assays -= 1
    if n_plate_blanks:
        assays_per_plate -= n_plate_blanks
    n_plates = int(np.ceil(n_assays / assays_per_plate))
    n_controls = n_plates * n_plate_controls if n_plate_controls else 0
    n_blanks = n_plates * n_plate_blanks if n_plate_blanks else 0
    assay_counts = utils.partition(n_assays + n_controls + n_blanks, n_plates)
    upto = lambda k: [i + 1 for i in range(k)]
    well_locs = [upto(ac * replicates) for ac in assay_counts]
    for locs in well_locs:
        locs.reverse()
        if randomize:
            random.shuffle(locs)

    control_locations = []
    if n_plate_controls:
        for i, locs in enumerate(well_locs):
            control_locations.append(
                (i + 1, [locs.pop() for _ in range(replicates * n_plate_controls)])
            )
    blank_locations = []
    if n_plate_blanks:
        for i, locs in enumerate(well_locs):
            blank_locations.append(
                (i + 1, [locs.pop() for _ in range(replicates * n_plate_blanks)])
            )

    assay_locations = []
    current = 0
    while current < n_assays:
        for i, locs in enumerate(well_locs):
            if current >= n_assays:
                break
            loc = (i + 1, [locs.pop() for _ in range(replicates)])
            assay_locations.append(loc)
            current += 1

    empty_locations = []
    for i, ac in enumerate(assay_counts):
        empty_locations.append(
            (i + 1, [k + 1 for k in range(ac * replicates, wells_per_plate)])
        )

    return dict(
        n_plates=n_plates,
        controls=control_locations,
        assays=assay_locations,
        blanks=blank_locations,
        empty=empty_locations,
    )


def layout_plates(
    solutions,
    wells_per_plate=96,
    replicates=1,
    plate_control=None,
    plate_blank=None,
    randomize=True,
):
    """Arrange replicates, controls, and blanks on multiwell plates.

    Inputs
    ------
    solutions: list of ingredients.Solution
        Solutions to be arrayed in plates.
    wells_per_plate: int, default=96
        Total number of wells available on each plate.
    replicates: int, default=1
        Number of times each assay will be assigned a location. Replicates for
        a single assay appear on the sample plate. If `randomize` is False,
        the replicates appear in adjacent wells; otherwise, the replicates are
        scattered across the plate.
    plate_control: tuple(str, int), default=None
        If the ID of one of the assays is given along with an int N, 
        it will be repeated on every plate as a control N times.
    plate_blank: tuple(str, int), default=None
        If the ID of one of the assays is given along with an int N, 
        it will be repeated on every plate as a blank N times.
    randomize: bool, default=True
        If True, the assays and replicates are assigned to random wells on the
        plate, although each assay's replicates remain on the same plate.

    Returns
    -------
    A pd.DataFrame with columns:
        plate: int
            Plate IDs from LoadPlate().id.
        well: int
            1-indexed array of well IDs.
        solution_id: str
            IDs from solution.id.
        replicates: int
            Replicate number (1-indexed).
        plate_control: bool
            True if the row is a control assay.
        plate_blank: bool
            True if the row is a blank assay.
    """
    solution_ids = [solution.id for solution in solutions]
    n_solutions = len(solution_ids)
    layout = pd.DataFrame(
        dict(
            plate=0,
            well=0,
            solution_id=solution_ids,
            replicate=1,
            plate_control=False,
            plate_blank=False,
        )
    )

    _INFO("Assigning %i replicate(s) per assay.", replicates)

    layout = [layout.copy() for _ in range(replicates)]
    for i in range(replicates):
        layout[i]["replicate"] = i + 1
    layout = pd.concat(layout, ignore_index=True)

    # assign assays to wells
    _INFO("Identifying well locations.")
    locs = _make_locations(
        n_solutions,
        wells_per_plate,
        replicates,
        n_plate_controls=plate_control[1] if plate_control else 0,
        n_plate_blanks=plate_blank[1] if plate_blank else 0,
        randomize=randomize,
    )

    _INFO("Final configuration has %i plates.", locs["n_plates"])

    if plate_blank is not None:
        _INFO("Using plate blanks.")
        blanks = layout[layout.solution_id == plate_blank[0]]
        accum = pd.DataFrame()

        for i in range(locs["n_plates"]):
            new = blanks.copy()
            new.plate_blank = True
            new.plate = locs["blanks"][i][0]
            new.replicate = list(range(1, replicates * plate_blank[1] + 1))
            keys = pd.Series(locs["blanks"][i][1])
            new.well = keys[new.replicate - 1].values
            accum = pd.concat([accum, new], ignore_index=True)
        blanks = accum
    else:
        blanks = pd.DataFrame()

    if plate_control is not None:
        _INFO("Using plate controls.")
        layout.plate_control = layout.solution_id == plate_control[0]
        controls = layout[layout.plate_control]
        others = layout[~layout.plate_control]
        accum = pd.DataFrame()
        for i in range(locs["n_plates"]):
            new = controls.copy()
            new.plate = locs["controls"][i][0]
            new.replicate = list(range(1, replicates * plate_control[1] + 1))
            keys = pd.Series(locs["controls"][i][1])
            new.well = keys[new.replicate - 1].values
            accum = pd.concat([accum, new], ignore_index=True)

        controls = accum
    else:
        controls = pd.DataFrame()
        others = layout

    _INFO("Assigning solutions to plates and wells.")
    others_idx = others.reset_index()
    others_grouped = others_idx.groupby(
        by=["solution_id", "plate_control", "plate_blank"],
    )

    others_aggregated = others_grouped.agg(
        replicate=("replicate", lambda x: len(x)), og_indexes=("index", tuple)
    )
    for _, row in others_aggregated.iterrows():
        location = locs["assays"].pop(0)
        for r in range(row["replicate"]):
            idx = row["og_indexes"][r]
            others.loc[idx, "plate"] = location[0]
            others.loc[idx, "well"] = location[1][r]

    layout = pd.concat([controls, blanks, others], ignore_index=True)

    return layout


# ===========================================================================
# =============================== scheduling ================================
# ===========================================================================


def _make_mantis_instructions(layout, volume, max_stocks, plate_type):
    stocks = set()
    for vols in volume.volume_added:
        stocks |= set(vols.index)
    stock_sets = utils.distribute(list(stocks), max_stocks)
    instructions = []
    for stocks in stock_sets:
        plate_instructions = []
        for plate in layout.plate.unique():
            pl = layout[layout.plate == plate]
            worklist = []
            for well, solution_id in zip(pl.well, pl.solution_id):
                vols = volume.loc[solution_id, "volume_added"]
                vols = vols[vols.index.intersection(stocks)]
                if not len(vols):
                    continue
                worklist.append(
                    pd.DataFrame(
                        dict(well=well, stock=vols.index, volume=[x for x in vols])
                    )
                )
            if not worklist:
                continue
            worklist = pd.concat(worklist)
            plate_instructions.append(
                (
                    LoadPlate(str(plate), plate_type.id_),
                    Pipette(worklist, makeids.unique_id_from_time(length=8)),
                )
            )
        instructions.append(
            RunMantis(stock=LoadStocks(sorted(stocks)), plates=plate_instructions)
        )

    return layout.plate.unique(), instructions


def schedule_mantis(
    solutions,
    strains,
    environments,
    stocks,
    plate=constants.labware["WP384"],
    replicates=1,
    plate_control=None,
    plate_blank=None,
    randomize=True,
    plate_prefix="",
    total_volume=None,
    working_volume=None,
    min_drop_size="0.1 ul",
    quantity_units="ug/ul",
    excess=None,
    max_stocks=24,
):

    solutions = utils.assert_list(solutions)
    strains_dict = dict()
    env_dict = dict()
    if isinstance(strains, dict):
        strains_dict = strains
        strains = strains.keys()
    strains = utils.assert_list(strains)
    if isinstance(environments, dict):
        env_dict = environments
        environments = environments.keys()
    environments = utils.assert_list(environments)

    if working_volume is None:
        working_volume = plate.working_volume
    if total_volume is None:
        total_volume = working_volume

    plate_offset = 0
    layout = pd.DataFrame()
    if strains_dict or env_dict:
        solutions_original = copy.deepcopy(solutions)
    for environment in environments:
        for strain in strains:
            if strains_dict:
                solutions_strain = strains_dict[strain]
                solutions_subset = [
                    s
                    for s in solutions
                    if (
                        s.id in solutions_strain
                        or s.id == plate_control[0]
                        or s.id == plate_blank[0]
                    )
                ]
            if env_dict:
                solutions_env = env_dict[environment]
                solutions = [
                    s
                    for s in solutions
                    if (
                        s.id in solutions_env
                        or s.id == plate_control[0]
                        or s.id == plate_blank[0]
                    )
                ]

            sub_layout = layout_plates(
                solutions,
                wells_per_plate=plate.n_wells(),
                replicates=replicates,
                plate_control=plate_control,
                plate_blank=plate_blank,
                randomize=randomize,
            )
            sub_layout.plate += plate_offset
            plate_offset = max(sub_layout.plate)
            sub_layout["strain"] = strain
            sub_layout["environment"] = environment
            layout = pd.concat([layout, sub_layout], ignore_index=True)
            if strains_dict or env_dict:
                solutions = copy.deepcopy(solutions_original)

    plate_names = [
        makeids.unique_id(prefix="Plate " + plate_prefix + str(i + 1))
        for i in range(max(layout.plate))
    ]
    layout.plate = [plate_names[i - 1] for i in layout.plate]
    unique_solution_ids = set()
    unique_solutions = list()
    for s in solutions:
        if s.id not in unique_solution_ids:
            unique_solution_ids.add(s.id)
            unique_solutions.append(s)
    volume = plan_mantis_pipetting(
        unique_solutions,
        stocks,
        total_volume=total_volume,
        working_volume=working_volume,
        excess=excess,
        min_drop_size=min_drop_size,
        quantity_units=quantity_units,
    )

    plates, mantis_instructions = _make_mantis_instructions(
        layout, volume, max_stocks, plate
    )
    return plates, mantis_instructions, layout
