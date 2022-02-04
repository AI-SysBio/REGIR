## Interacting with SBML files

The Systems Biology Markup Language (SBML) is a representation format based on XML, used to store and share computational models. It is a standard for representing computational models in systems biology. REGIR can read and write models in the SBML format. If you want to wotk wit hSBML files, you first need to install the `libSBML` and `simplesbml` libraries with:

	- pip install libSBML
	- pip install simplesbml
	
	
Then, you can store any of your models builtwith REGIR with the `get_model_in_SBML()` method in REGIR. In REGIR, the rate laws are always proportional to the 

    
Then, you can store any of your models builtwith REGIR with the `build_model_from_SBML` method in REGIR:

class param:
    Tend = 200		#Length of the simulation
    unit = 'h'		#Unit of time (is used for plotting purpose only)
    N_simulations = 20	#The simulation results should be averaged over many trials
    timepoints = 100	#Number of timepoints to record (make surethat this number isnt too big)
    


def main():
    
    #https://www.ebi.ac.uk/biomodels/BIOMD0000000478#Files
    SBML_file = 'stem_cell_model.sbml'
    REGIR_SBML = REGIR_from_SBML(verbose=True)
    G_simul = REGIR_SBML.build_model_from_SBML(SBML_file, param)
 
While many SBML models are available in the BioModels database (https://www.ebi.ac.uk/biomodels/), the majority of them are unfortunatly not compatible with REGIR. In fact, REGIR currently only supports rate laws that are proportional to reactant amounts and that do not depend on other populations. Also, As REGIR 
