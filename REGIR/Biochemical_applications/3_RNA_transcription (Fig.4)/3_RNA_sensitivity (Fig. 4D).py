import numpy as np
import matplotlib.pyplot as plt
from SALib.sample import saltelli
from SALib.analyze import sobol 
#https://salib.readthedocs.io/en/latest/

from RNA_model import Gillespie_kinetics



N_simulations = 200
sensitivity_res = 64
recompute_sensitivity = False


class param_opt:
    n_iteration = 100
    current_iter = 0
    use_already_computed_results = True
    
    
def main():
    
    #2- sensitivity analysis
    r_on = 1
    r_off = 1
    r_prod = 50

    Sobol_sensitivity_analysis(r_on, r_off, r_prod)
        
        
def Sobol_sensitivity_analysis(r_on, r_off, r_prod):
    
    """
    sensitivity of the coefficient of variation (CV) of mature mRNA 
    (SD divided by the mean) to small perturbations.
    The CV is averaged over the cell cycle.
    
    Sensitivity of a parameter r is defined formally as Sr = r/CV * dCV/dr,
    meaning that a 1% change in the value of parameter r leads to Sr% change in CV
    
    See https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis
    """
      
    
    rmin = min(1,np.min([r_on, r_off, r_prod]))
    Tend = 5/rmin

    problem = {
      'num_vars': 6, 
      'names': [r'$\lambda_{on}$', r'$\lambda_{off}$', r'$\lambda_{prod}$', r'$\alpha_{on}$', r'$\alpha_{off}$', r'$\alpha_{prod}$'], 
      'bounds': [[r_on*0.7, r_on*1.3],[r_off*0.7, r_off*1.3],[r_prod*0.7, r_prod*1.3],[2, 10],[2, 10],[2, 10]]
    }
    
    if recompute_sensitivity:
        # Generate samples
        param_values = saltelli.sample(problem, sensitivity_res, calc_second_order=False)
        n_sim_tot = len(param_values)
        
        def evaluate_Gillespie(values):
            Y = np.empty((values.shape[0],4))
            for i, X in enumerate(values):
                print('  Simulation %s/%s ...' % (i+1,n_sim_tot))
                param_list = dict()
                param_list['r_deg'] = 1
                param_list['r_on'] = X[0]
                param_list['r_off'] = X[1]
                param_list['r_prod'] = X[2]
                param_list['alpha_on'] = X[3]
                param_list['alpha_off'] = X[4]
                param_list['alpha_prod'] = X[5]
                param_list['Tend'] = Tend
                print('    r_on: %.2f, r_off: %.2f, r_prod: %.2f, \n    alpha_on: %.2f, alpha_off: %.2f, alpha_prod: %.2f' % (param_list['r_on'],param_list['r_off'],param_list['r_prod'],param_list['alpha_on'],param_list['alpha_off'],param_list['alpha_prod']))
                
                mean, CV, Entropy, Fano, _ = Gillespie_kinetics(param_list, plot = False, N_simulations = N_simulations)
                Y[i,:] = np.array([mean, CV, Entropy, Fano])
                
            return Y
        
        # Run model on all values
        Y = evaluate_Gillespie(param_values)
        
        np.save('sensitivity_analysis.npy',Y)
        
    else:
        Y = np.load('sensitivity_analysis.npy')
    
    measures = ['mean','CV','Fano','Entropy']
    for i in range(4):
        # Perform analysis
        Si = sobol.analyze(problem, Y[:,i], calc_second_order=False)
        # Returns a dictionary with keys 'S1', 'S1_conf', 'ST', and 'ST_conf'
        # (first and total-order indices with bootstrap confidence intervals)
        
        print('  %s sensitivity analysis:' % measures[i])
        print('   ', Si['S1'])
        print('   ', Si['S1_conf']) 
        
        plt.figure(figsize = (5.5,4))
        plt.errorbar(np.arange(problem['num_vars']), Si['S1'], yerr=Si['S1_conf'], fmt =  'o', color = 'black', capsize=5, elinewidth=2, markeredgewidth=2, markersize = 10, markerfacecolor = 'red')
        plt.xticks(ticks = np.arange(problem['num_vars']), labels = problem['names'], rotation = 90, size = 19)
        plt.axhline(y=0, color = 'black', linestyle = '--')
        plt.xlim(-0.5,5.5)
        plt.ylabel('Variance-based Sensitivity')
        plt.show()
    

    
if __name__ == "__main__":
    plt.rcParams.update({'font.size': 16})
    main()