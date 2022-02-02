"""
  install:
    pip install simplesbml
    pip install tellurium
    pip install REGIR
    conda install libSBML
    
  See https://simplesbml.readthedocs.io/en/latest/
"""



from REGIR_SBML import REGIR_from_SBML


#Set the simulation parameters:
class param:
    Tend = 200		#Length of the simulation
    unit = 'h'		#Unit of time (is used for plotting purpose only)
    N_simulations = 20	#The simulation results should be averaged over many trials
    timepoints = 100	#Number of timepoints to record (make surethat this number isnt too big)
    


def main():
    
    """
     Shape_and_distrib_in_local_param:
        
        For non markovian SBML models, the shape parameters are stored in the local paramater of each reaction.
        However, if your SBML file is from another source (such as BioModel), there will be no information about 
        the shape parameter and thus you will need to define it yourself after loading the SBML file
    """
    
    
    SBML_file = 'stem_cell_model.sbml'
    G_simul = REGIR_from_SBML(SBML_file, param, Shape_and_distrib_in_local_param = True)
    
    
    print('Running %s REGIR simulation for Tend = %s%s ...' % (param.N_simulations,param.Tend,param.unit))    
    G_simul.run_simulations(param.Tend, verbose = False)
    #G_simul.plot_inter_event_time_distribution()
    G_simul.plot_populations()
    

    
    
if __name__ == '__main__':
    main()
