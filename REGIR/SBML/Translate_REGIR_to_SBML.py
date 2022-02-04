
import REGIR as gil


""" --------------------------------------------------------------------------------------------
2 reactions:
    ESC -> EPI,  differentiation
    EPI -> NPC,  differentiation
"""

class param:
    Tend = 170
    unit = 'h'
    N_simulations = 20
    timepoints = 100
    

def main():
    
    print("\n=================================================================================") 
    print("=========================== Rejection Gillepsie algorithm =========================") 
    print("=================================================================================\n\n")                 

    r_diffAB = 0.0339
    r_diffBC  = 0.0164
    alpha_diffAB = 27.60
    alpha_diffBC =  39.29
    

    #initialise reactants
    N_init = dict()
    N_init['ESC'] = 100
    N_init['EPI'] = 0
    N_init['NPC'] = 0
    
    #initialise reaction channels
    EE_differentiation = gil.Reaction_channel(param,rate=r_diffAB, shape_param=alpha_diffAB, distribution = 'Gamma', name='Differentiation: ESC -> EPI')
    EE_differentiation.reactants = ['ESC']
    EE_differentiation.products = ['EPI']
    
    EN_differentiation = gil.Reaction_channel(param,rate=r_diffBC, shape_param=alpha_diffBC, distribution = 'Gamma', name='Differentiation: EPI -> NPC')
    EN_differentiation.reactants = ['EPI']
    EN_differentiation.products = ['NPC']
    reaction_channel_list = [EE_differentiation,EN_differentiation]
    
    #initialise the Gillespie simulation
    G_simul = gil.Gillespie_simulation(N_init,param)
    G_simul.reaction_channel_list = reaction_channel_list
    print(G_simul)
    
    sbml_model = G_simul.get_model_in_SBML()
    print(sbml_model.toSBML())
    with open("stem_cell_model.sbml", "w") as text_file:
        text_file.write(sbml_model.toSBML())
    
    
if __name__ == '__main__':
    main()