
import REGIR_test as gil
import sys
import numpy as np
import simplesbml

"""
   install:
     pip install simplesbml
     pip install tellurium
     pip install REGIR
     pip/conda install libSBML
        
     See https://simplesbml.readthedocs.io/en/latest/
"""

class REGIR_from_SBML:
    
    def __init__(self, verbose = False):
        self.verbose = verbose

            
    def build_model_from_SBML(self,SBML_file, param, Shape_and_distrib_in_local_param = False):


    
        print('Extraction of SBML model from %s...' % SBML_file)
        print()
        
    
    
        self.model = simplesbml.loadSBMLFile(SBML_file)
        
        #1 = Get relevant parameters from the SBML model
        self.fspecies = self.model.getListOfFloatingSpecies()
        self.bspecies = self.model.getListOfBoundarySpecies()
        self.species = self.fspecies + self.bspecies
        self.reactions = self.model.getListOfReactionIds()
        self.compartements = self.model.getListOfCompartmentIds()
        self.functions = self.model.getListOfFunctionIds()
     
        """         
        print ('List of rules = ', Function_IDs)
        for func in functions:
            print(model.getFunctionBody(func))
            print(model.getNumArgumentsInUserFunction(func))
            print(model.getListOfArgumentsInUserFunction(func))
            
      
        Rule_IDs = model.getListOfRuleIds()
        print ('List of rules = ',Rule_IDs)
        for rule in Rule_IDs:
            print(model.getRuleRightSide(rule))
            print(model.getRuleType(rule))
        """
       
            
        #sys.exit()
        
        #2 Incorporate into the REGIR framework 
        N_init = dict()
        
        
        Concentration_array = np.zeros(len(self.species))
        for ei,entity in enumerate(self.species):
            if self.model.isAmount(entity):
                n_init = int(self.model.getSpeciesInitialAmount(entity))
                N_init[entity] = n_init
                if n_init > 10000:
                    print(" Warning: REGIR be very slow if initial amount is very high")
                if n_init > 100000 and False:
                    print(" Simulating more than 100k entities with REGIR is unrealistic")
                    sys.exit()
            elif self.model.isConcentration(entity):
                concentration = self.model.getSpeciesInitialConcentration(entity)
                Concentration_array[ei] = concentration
            else:
                #print('Entity %s is neither given in Amount or concentration, we set 0' % entity)
                N_init[entity] = 0
                
        if np.sum(Concentration_array) > 0:
            Concentration_array = Concentration_array/np.min(Concentration_array[Concentration_array>0])*10
            for ei,entity in enumerate(self.species):
                if self.model.isConcentration(entity):
                    N_init[entity] = int(Concentration_array[ei])
    
        if self.verbose: print('  Initial population:',N_init)
        
            
        reaction_channel_list = []
        self.parameter_list = sorted(self.model.getListOfParameterIds(), key=len)[::-1]
        for ri,reactionId in enumerate(self.reactions):
            
            if self.verbose: print()
            
            reactants_list = []
            n_reactants = self.model.getNumReactants(reactionId)
            for nri in range(n_reactants):
                reactants_list.append(self.model.getReactant(reactionId, nri))
            products_list = []
            n_products = self.model.getNumProducts(reactionId)
            for pi in range(n_products):
                products_list.append(self.model.getProduct(reactionId, pi))         
            if self.verbose: print("  Reaction", reactants_list, "->", products_list,":")
            
            
            rate_law = self.model.getRateLaw(reactionId)
            Reactions_to_create, rate = self.analyze_and_transform_rate_law(rate_law)
            
            
            
            if Shape_and_distrib_in_local_param:
                local_param = get_local_param(self.model,ri)
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
        
        if self.verbose: print()
        print("Extracted", G_simul)
    
        
        return G_simul


    def analyze_and_transform_rate_law(self,rate_law):
        
        
        #First, replace the function by its actual expression
        for func in self.functions:
            if func in rate_law:
                if self.verbose: print('  ',func)
                func_arguments = self.model.getListOfArgumentsInUserFunction(func)
                func_body = self.model.getFunctionBody(func)
                rate_law = evaluate_functions(rate_law,func,func_body,func_arguments)
        
        if self.verbose: print('  ',rate_law)   
        
        #Then replace by the parameters
        for paramID in self.parameter_list:
            rate_law = rate_law.replace(paramID,str(self.model.getParameterValue(paramID)))
            
        for ci,comp in enumerate(self.compartements):
            comp_size = get_comp_size(self.model,ci)
            comp_size = 1 #Let's say compsize is always 1
            rate_law = rate_law.replace(comp,str(comp_size))
            
            
        if self.verbose: print('  ',rate_law)   
        
        #Detect what reactions should be created (+,-,/,*, or incompatible)
        
            
        for entity in self.species:  #For REGIR implementation, the rate is always proportional to the number of reactants
            rate_law = rate_law.replace(entity,'1')

        
        try:
            rate = eval(rate_law)
        except Exception as e:
            print('  ERROR: %s' %e)
            print('  When computing rate = ',rate_law)
            sys.exit()
        if self.verbose: print("   rate = %s" % rate)
    
    
    
        reactions_to_create = None
        return reactions_to_create, rate


def evaluate_functions(rate_law,func,func_body,func_arguments):
    
    full_func = func + rate_law.split(func)[1].split(')')[0]+')'
    func_args = rate_law.split(func)[1].split(')')[0].replace('(','').replace(' ','').split(',')
    for ai,arg in enumerate(func_args):
        func_body = func_body.replace(func_arguments[ai],func_args[ai])
    rate_law = rate_law.replace(full_func,func_body)
    
    return rate_law


def get_local_param(model,ri):
    v = model.getModel().getListOfReactions()[ri]
    local_ids = []
    local_values = []
    for k in v.getKineticLaw().getListOfParameters():
        local_ids.append(k.getId())
        local_values.append(k.getValue())
    local_params = dict(zip(local_ids, local_values))
    return local_params


def get_comp_size(model,ci):
    comp = model.getModel().getListOfCompartments()[ci]
    return comp.getSize()


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
        