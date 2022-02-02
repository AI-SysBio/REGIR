""" 
==============================================================================================================
================================= Below is only REGIR code related to SBML ===================================
==============================================================================================================
"""

import REGIR_test as gil
import sys
import simplesbml

"""
   install:
     pip install simplesbml
     pip install tellurium
     pip install REGIR
     pip/conda install libSBML
        
     See https://simplesbml.readthedocs.io/en/latest/
"""

def REGIR_from_SBML(SBML_file, param, Shape_and_distrib_in_local_param = False):
    
    print('Extraction of SBML model from %s...' % SBML_file)
    print()
    


    model = simplesbml.loadSBMLFile(SBML_file)
    
    #1 = Get relevant parameters from the SBML model
    species = model.getListOfFloatingSpecies()
    reactions = model.getListOfReactionIds()
    
    #2 Incorporate into the REGIR framework 
    N_init = dict()
    N_av = 6.02214076e23
    for entity in species:
        if model.isAmount(entity):
            n_init = int(model.getSpeciesInitialAmount(entity))
        elif model.isConcentration(entity):
            n_init = int(model.getSpeciesInitialConcentration(entity)*N_av)
            print(" Warning: The Gillespie algorithmcan be very slow if initial amount is very high")
        else:
            sys.exit()
        N_init[entity] = n_init
        
    #print(N_init)
        
    reaction_channel_list = []
    parameter_list = model.getListOfParameterIds()
    for ri,reactionId in enumerate(reactions):
        
        reactants_list = []
        rate_law = model.getRateLaw(reactionId)
        n_reactants = model.getNumReactants(reactionId)
        for nri in range(n_reactants):
            reactants_list.append(model.getReactant(reactionId, nri))
            rate_law = rate_law.replace(model.getReactant(reactionId, nri),'1')
        products_list = []
        n_products = model.getNumProducts(reactionId)
        for pi in range(n_products):
            products_list.append(model.getProduct(reactionId, pi))         
        #print("reaction", reactants_list, "->", products_list,":")

        
        for paramID in parameter_list:
            rate_law = rate_law.replace(paramID,str(model.getParameterValue(paramID)))
        
        rate = eval(rate_law)
        #print("  rate = %s" % rate)
        
        if Shape_and_distrib_in_local_param:
            local_param = get_local_param(model,ri)
            shape_param = local_param['shape_param']
            distribution_index = int(local_param['distribution_index'])
            distribution = get_distribution(distribution_index)
        else:
            shape_param = 1
            distribution = 'exp'

        
        reaction = gil.Reaction_channel(param,rate=rate, shape_param=shape_param, distribution = distribution, reactants = reactants_list, products = products_list, name = reactionId)
        reaction_channel_list.append(reaction) 
        
        
        #local_param = model.getlocal_param(reactionId)
        #print(local_param)
        

        
    G_simul = gil.Gillespie_simulation(N_init,param)  
    G_simul.reaction_channel_list = reaction_channel_list
    
    print("Extracted", G_simul)

    
    return G_simul


def get_local_param(model,ri):
    
    v = model.getModel().getListOfReactions()[ri]
    local_ids = []
    local_values = []
    for k in v.getKineticLaw().getListOfParameters():
        local_ids.append(k.getId())
        local_values.append(k.getValue())
    local_params = dict(zip(local_ids, local_values))
    
    return local_params
    



def get_distribution(index):
    distribution_list = ['exp','weib','norm','gamma','lognorm','cauchy']
    return distribution_list[index]

def Distribution_index(distribution):
    if distribution.lower() in ['exponential', 'exp']:
        return 0                    
    elif distribution.lower() in ['weibull','weib']:
        return 1                            
    elif distribution.lower() in ['gaussian', 'normal', 'norm']:
        return 2                           
    elif distribution.lower() in ['gamma','gam']:
        return 3                          
    elif distribution.lower() in ['lognormal','lognorm']:
        return 4                       
    elif distribution.lower() in ['cauchy', 'cau']:
        return 5
        