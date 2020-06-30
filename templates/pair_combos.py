"""pair_combos.py"""

import datetime 
import os

from plateplan import (
	constants, ingredients, makeids, mantis, mapper, scheduling, units
)

amino_acid_ids = ["ALA", "ARG", "ASP"]
amino_acid_names = ["alanine", "arginine", "aspartic acid"]
final_conc = {
	"ALA": "10 mg/ml",
	"ARG": "10 mg/ml",
	"ASP": "10 mg/ml"
}

reagents = {}
stocks = []

# Define stock solutions and reagents
for i in range(3):
	id_ = amino_acid_ids[i] # Get id value
	name = amino_acid_names[i] # Get name value
	conc = str(units.parse(final_conc[id_]) * 100) #Get 100x stock conc value

	reagents[id_] = ingredients.Reagent(id_=id_,  name=name) # Not used in this experiment
	soln = ingredients.Solution(reagents={id_: conc},  id_=id_+"_stock") # Make solution
	print(soln)
	stocks.append(ingredients.Stock(
		ingredient=soln,
		date_made=datetime.date.today(), # Today's date
		date_expires=datetime.date.today() + datetime.timedelta(6  *  365  /  12), # +6 months for expiry date
		quantity=units.parse("10 ml"), # Amount of stock we have stored
		labware="Conical15ml", # What its stored in
	)) # Make stock from Solution, soln and add it to stocks list

solutions = []
# Define solutions you want to make
pairs = [["ALA", "ARG"], ["ALA", "ASP"], ["ARG", "ASP"]]
for first, second in pairs:
	reagents = {
		first: final_conc[first],
		second: final_conc[second]
	}
	id_ = first + "_" + second + "_mixture"
	soln = ingredients.Solution(reagents=reagents,  id_=id_) # Make solution
	solutions.append(soln)

# Add water to solutions so we can use it as a blank
water = ingredients.Solution(reagents={"water": "0 g/l"},  id_="water") # Make solution
solutions.append(water)

# Schedule the experiments with your parameters
plates, instructions, layout = scheduling.schedule_mantis(
    solutions,
    [constants.strains["SMU"]],
    ["aerobic"],
    stocks,
    excess="water",
    plate=constants.labware["WP384"],
    total_volume="80 ul",
    working_volume="78 ul",
    plate_control=("ALA_ARG_mixture",  5),
    plate_blank=("water",  5),
    replicates=10, 
    min_drop_size="0.1 ul"
)


# Make an InstructionSet based on the instructions outputted from
# the scheduler above. Used for generating the files below.
friendly_name = "Pairs"
exp_id = makeids.unique_id(prefix=friendly_name)
instruction_set = scheduling.InstructionSet(
	id_=exp_id, instructions=instructions, owner="Adam",
)

# Export the files needed for the experiment
output_filepath = os.path.join("experiment_outputs", exp_id)
os.makedirs(output_filepath)

# Export files in Mantis-friendly format
mantis.generate_experiment_files(
	instruction_set, layout, "2 ul", path_to_worklists=output_filepath # 2 ul dispenses of the bacteria strain per well 
)
# Export an experiment map for future use when mapping data back to each experiment
mapper.save_well_map(
	path=output_filepath,
	layout=layout,
	instruction_set_id=instruction_set.id,
	dimensions=constants.labware["WP384"].shape,
) # Saves map.csv