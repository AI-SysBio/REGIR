
from REGIR_SBML import REGIR_from_SBML


#Set the simulation parameters:
class params:
    Tend = 200		    #Length of the simulation
    N_simulations = 20
    unit = 'h'		    #Unit of time (used for plotting only)
    timepoints = 100	#Number of timepoints to record (used for plotting only))


def main():
   
    SBML_file = 'stem_cell_model.sbml'
    REGIR_SBML = REGIR_from_SBML(verbose=True)
    G_simul = REGIR_SBML.build_model_from_SBML(SBML_file, params)
    
    print('Running %s REGIR simulation for Tend = %s%s ...' % (params.N_simulations,params.Tend,params.unit))    
    G_simul.run_simulations(params.Tend, verbose = False)
    #G_simul.plot_inter_event_time_distribution()
    G_simul.plot_populations()
    

if __name__ == '__main__':
    main()
