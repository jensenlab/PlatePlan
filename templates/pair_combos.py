"""pair_combos.py"""

from jendrop import constants, ingredients, scheduler, units

amino_acid_ids = ["ALA", "ARG", "ASP"]
amino_acid_names = ["alanine", "arginine", "aspartic acid"]
final_conc = {
	"ALA": "10 mg/ml",
	"ARG": "10 mg/ml",
	"ASP": "10 mg/ml"
}

reagents = {}
stocks = {}

# Define stock solutions and reagents
for i in range(3):
    id_ = amino_acid_ids[i] #Get id value
    name = amino_acid_names[i] #Get name value
    conc = str(units.parse(final_conc[id_]) * 100) #Get 100x stock conc value
    
	reagents[id_] = ingredients.Reagent(id_=id_,  name=name) #Not used in this experiment
	soln = ingredients.Solution(reagents={id_: conc},  id_=id_+"_stock") #Make solution
	
	stocks[id_] =  ingredients.Stock(
		ingredient=solution,
		date_made=datetime.date.today(), #Today's date
		date_expires=datetime.date.today() + datetime.timedelta(6  *  365  /  12), #+6 months for expiry date
		quantity=units.parse("10 ml"), #Amount of stock we have stored
		labware="Conical15ml", #What its stored in
	) #Make stock from Solution, soln

solutions = []

# Define solutions you want to make
pairs = [["ALA", "ARG"], ["ALA", "ASP"], ["ARG", "ASP"]]
for first, second in pairs:
	reagents = {
		first: final_conc[first],
		second: final_conc[second]
	}
	id_ = first + "_" + second + "_mixture"
	soln = ingredients.Solution(reagents=reagents,  id_=id_) #Make solution
	solutions.append(soln)

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
    plate_control=("CDM",  1),
    plate_blank=("CDM",  1),
    replicates=4, 
    min_drop_size="0.1 ul"
)