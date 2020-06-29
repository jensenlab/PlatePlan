
![PlatePlan Logo](img/plateplan_logo.png)
# Welcome to PlatePlan
PlatePlan is meant to be the bridge between your hypotheses and results. In order to make that jump, a lot goes into the logistics and planning out an experiment, from how many multi-well plates to use and where to dispense each reagent, to which reagents to batch together and how much volume to dispense. This logistics are especially challenging in the case of high-throughput assays where you may have tens of thousands of experiments to manage and keep track of. 

PlatePlan takes some of that burden away from the experimenter by handling those logistics for you. It has been tailored to our lab's set up ([Jensen Lab, UIUC](www.jensenlab.net)) to take advantage of the [Formulatrix Mantis](https://formulatrix.com/liquid-handling-systems/mantis-liquid-handler/) liquid handling robot for combinatorial bacteria growth assays, but can be adapted to fit another use case with relative ease. 

It does so by generalizing an experimental structure is explained further in depth in **[Structure](#structure)**. That being said, we've tried to make it easy for an experimenter to get started. To make an experiment, the biggest hurdle to use PlatePlan is determining what final solutions you want your liquid handler to complete for you--this comes in the form of a list of Solutions that each map a Reagent to its final concentration. Beyond inputting a few extra experimental parameters like how many replicates or blanks you need, that's it.

You can refer to **[Running an Experiment](#running-an-experiment)** for more detailed requirements and instructions or look at our sample experiment in `templates/l2o_amino_acids.py` which is an experiment we've developed to perform a large scale leave-two-out growth assay for amino acids in a chemically defined medium ([CDM](https://en.wikipedia.org/wiki/Chemically_defined_medium)). A simpler example experiment is also located in `templates/pair_combos.py` and used as the example code for the sections below.

# Software requirements
1. A copy of [Gurobi optimizer](https://www.gurobi.com/) must be install in your environment, with the [Python interface](https://www.gurobi.com/documentation/quickstart.html) set up.
2. Clone this repo using `git clone [INSERT REPO LINK HERE]`
3. Pipenv environment for this project
	- [install pipenv](https://pipenv-fork.readthedocs.io/en/latest/install.html)
	- run `pipenv install` in the root project directory to install the required dependencies
4. 

# Structure
In order to define an experiment, the user needs to define the input Solutions, Reagents, and Stocks. Stocks can either be of Solution or Reagent type--they are the physical analog of Solutions and Reagents which are purely virtual themselves. 

## Ingredients
### Reagents
Reagents are the smallest unit of the experiment you are planning. They have two 
properties: 
`id`: An identifier string of the reagent, e.g. 'sodium_chloride'
`name`: A secondary name of the reagent, e.g. 'NaCl'
`molecular_weight`: A string of the molecular weight, e.g. '58.44 g/mol'
### Solutions
Solutions are composed of Reagent.id_ and their concentrations in solution and represent both the solutions you need to dispense and the solutions you are trying to make. They are composed of 3 components:
properties: 
`id_`: A string of the solution name, e.g. '2x sodium chloride'
`reagents`: A dictionary mapping Reagent.name -> concentration (string), e.g. 

       {
	       'sodium_chloride' : '50 mg/L',
	       'sodium_carbonate' : '50 mg/L'
	   }

`solvent`: A string of the solvent

### Stocks
Stocks contain a single Ingredient (Reagent or Solution) with additional real-world laboratory information, useful for inventory and management. They are composed of the following attributes:
`ingredient`: The corresponding Reagent or Solution object
`id_`: The ID of this new stock
`quantity`: The volume of the stock you have stored
`labware`: The type of labware the stock is stored in
`location`: The location of the stock
`date_made`: The date the stock was made
`date_expires`: The expiry date of the stock

### Example experiment set up using stocks and solutions
For this example experiment, we will make all combinations of 2 of the follow three amino acids:
**alanine, arginine, aspartic acid**. We have them stored as 100x stocks in our lab. And the final concentration of each will be 10 mg/ml. We need to generate 3 Reagents (alanine, arginine, aspartic acid) , 3 Reagent Stocks (alanine stock, arginine stock, aspartic acid stock), and 3 Solutions to make (alanine+arginine, alanine+aspartic acid, arginine+aspartic acid).

#### Generate 3 Solution Stocks:
	"""pair_combos.py"""
	from plateplan import constants, ingredients, scheduler, units
	
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
	
	...
	
#### Generate Solutions to Make:
	"""pair_combos.py"""
	...
	
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
		
	...

# Running an Experiment
Once you have generated all the necessary Solutions and Stocks needed to run the experiment, you must now define the other experiment parameters and run the scheduler.
### Defining your experimental parameters
`solutions`: list of ingredients.Solutions you want to make
`strains`: list of str of the bacteria you want to use
`environments`: list of str of the environments you want to grow in
`stocks`: list of ingredients.Stocks that you have for dispensing
`plate`: a constants.labware Labware plate, e.g. constants.labware["WP384"] for 384-well plate
`replicates`: int number of replicates 
`plate_control`: tuple of (Solution.id_, quantity) to use for the control, e.g. ("CDM",  1)
`plate_blank`:  tuple of (Solution.id_, quantity) to use for the blank, e.g. ("CDM",  1)
`randomize`: boolean flag to set if replicates should have randomized locations within the plate (replicates will always be placed on the same plate)
`plate_prefix`: str to add to the plate name
`total_volume`: str of volume in each well, e.g. "80 ul"
`working_volume`: str of the working volume in each well, e.g. "80 ul"
`min_drop_size`: str of the volume corresponding to the maximum resolution for your liquid handler, e.g. "0.1 ul"
`quantity_units`: str of concentration unit to use for liquid handler, e.g. "ug/ul"
`excess`: str of reagent to use to input for excess volume/to get to the correct concentrations, e.g. "water"
`max_stocks`: int of the maximum number of stocks you can load at a time in your liquid handler, e.g. 24 (we use a Formulatrix Mantis + LC3)

### Example scheduler run
	"""pair_combos.py"""
	
	...
	
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
	
	...

# Tools
We have created a few tools for the user and/or the scheduler.
## Constants 
Makes conversions between units easier and automatic. It is used by internally by the scheduler and by the user when defining experiments. You are able to perform arithmetic and conversions between magnitudes of the same unit and into other units types entirely.
## Mantis 
Converts from the scheduler plate plans to worklists compatible with the Formulatrix Mantis. 
## Mapper
Maps raw data back to the experiments planned by the scheduler. To do this, you need the plate plan map.csv that you generated and the raw files. Currently there is only support for reading in CSV data from a BioTek Epoch plate reader for our use case, but you can easily create a compatibility layer of your own
>The map.csv contains plate # and well #. 
# Publications
Please reach out if you use this software in any of you questions at manager@jensenlab.net, we'd love to hear about your project.

# Usage Policy/License
This software is provided to you with no guarantees under the [INSERT LICENSE NAME HERE].