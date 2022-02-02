"""
  install:
    pip install simplesbml
    pip install tellurium
    conda install libSBML
    
  See https://simplesbml.readthedocs.io/en/latest/
"""
import simplesbml

from REGIR_SBML import Distribution_index


#Define params
r_ESC_EPI = 0.032
r_EPI_NPC  = 0.017
alpha_ESC_EPI = 27.60
alpha_EPI_NPC =  39.29
ESC_EPI_distribution = "gamma"
EPI_NPC_distribution = "gamma"
ESC_EPI_params = {'shape_param':alpha_ESC_EPI, 'distribution_index':Distribution_index(ESC_EPI_distribution)}
EPI_NPC_params = {'shape_param':alpha_EPI_NPC, 'distribution_index':Distribution_index(EPI_NPC_distribution)}



#Define model
model = simplesbml.SbmlModel()
model.addCompartment(1, comp_id='comp')
model.addSpecies('ESC', 100, comp='comp')
model.addSpecies('EPI', 0, comp='comp')
model.addSpecies('NPC', 0, comp='comp')
model.addParameter('r1', r_ESC_EPI)
model.addParameter('r2', r_EPI_NPC)
model.addReaction(['ESC'], ['EPI'], 'ESC*r1', local_params = ESC_EPI_params, rxn_id='ESC_diff')
model.addReaction(['EPI'], ['NPC'], 'EPI*r2', local_params = EPI_NPC_params, rxn_id='EPI_diff')

"""
  model.addReaction(reactants, products, expression, local_params={}, rxn_id='')

  -> reactants and products are lists of species ids that the user wishes to 
     define as reactants and products, respectively. 

  -> Expression is a string that represents the reaction rate expression.

  -> local_params is a dictionary where the keys are local parameter ids and the 
     values are the desired values of the respective parameters.
"""

#save model into file
print (model.toSBML())
with open("stem_cell_model.sbml", "w") as text_file:
    text_file.write(model.toSBML())



