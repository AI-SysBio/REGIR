import numpy as np
import pandas as pd
import random as rd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from math import log, gamma
from scipy.stats import expon, weibull_min
from sklearn.metrics import mean_squared_error
import os,sys

import REGIR as gil


""" --------------------------------------------------------------------------------------------
2 reactions:
    ESC -> EPI,  differentiation
    EPI -> NPC,  differentiation
    
"""

recompute_all = False
sigmoid_adjusted = False


class param:
    Tend = 200
    unit = 'h'
    rd.seed(101)                #length of the simulation, in hours (41 days)
    N_simulations = 100          #The simulation results should be averaged over many trials
    timepoints = 170            #Number of timepoints to record (make surethat this number isnt too big)
    distribution = 'Gamma' #'Lognormal", 'Gaussian', 'Gamma', 'Weibull'
    
    
class param_opt:
    n_iteration = 100
    current_iter = 0

def main_optimization():
    #Gillespie_RMMSE()
    
    distribution_list = ['Exponential','Lognormal', 'Gaussian', 'Gamma', 'Weibull']
    
    if recompute_all:
        for distr in distribution_list:
        
            param.distribution = distr
            
            param_opt.current_iter = 0
            Pop, rmse = LIPO_Optimization()
            Gillespie_RMMSE(*Pop)
            if sigmoid_adjusted:
                np.save('Optimization/Optimal_parameters_sig_%s.npy' % distr, Pop)
            else:
                np.save('Optimization/Optimal_parameters_%s.npy' % distr, Pop)
            print(Pop)
    
    for distr in distribution_list:
        
        print()
        print()
        print('  Optimal RMSE for %s:' % distr)
        if sigmoid_adjusted:
            Pop = np.load('Optimization/Optimal_parameters_sig_%s.npy' % distr)
        else:
            Pop = np.load('Optimization/Optimal_parameters_%s.npy' % distr)
        param.distribution = distr
        Gillespie_RMMSE(*Pop, plot = False)
        

def LIPO_Optimization():
    
    import dlib
        
    """
    Requires dlib python library, install with: 
        https://medium.com/analytics-vidhya/how-to-install-dlib-library-for-python-in-windows-10-57348ba1117f
    The optimization uses the maxLIPO algorithm:
        http://blog.dlib.net/2017/12/a-global-optimization-algorithm-worth.html
        
    experimental_data is a dict of dict that contains 
    population value at each time for each population:
        dict[reactant_str][time] = value
    """
    
    
    
    #parameters tooptimize and their bound
    param_bounds = dict()

    param_bounds['r_diffAB'] = [0.031 , 0.034]
    param_bounds['r_diffBC'] = [0.015 , 0.0175]
    
    if param.distribution.lower() in ['gamma']:
        param_bounds['alpha_diffAB'] = [1 , 50]
        param_bounds['alpha_diffBC'] = [1 , 50]
    elif param.distribution.lower() in ['weibull']:
        param_bounds['alpha_diffAB'] = [1 , 5]
        param_bounds['alpha_diffBC'] = [1 , 5]
    elif param.distribution.lower() in ['lognormal','gaussian','normal']:
        param_bounds['alpha_diffAB'] = [0.05 , 0.5]
        param_bounds['alpha_diffBC'] = [0.05 , 0.5]
    elif param.distribution.lower() in ['exponential']:
        function_opt = lambda r1,r2: Gillespie_RMMSE(r1,r2, plot = False, plot_RMSE = False)
        xopt,yopt = dlib.find_min_global(function_opt,
                           [param_bounds[key][0] for key in param_bounds],
                           [param_bounds[key][1] for key in param_bounds],
                           param_opt.n_iteration)
    
        return (xopt,yopt)
    else:
        sys.exit()
    
    
    function_opt = lambda r1,r2,alpha1,alpha2: Gillespie_RMMSE(r1,r2,alpha1,alpha2, plot = False, plot_RMSE = False)
    xopt,yopt = dlib.find_min_global(function_opt,
                           [param_bounds[key][0] for key in param_bounds],
                           [param_bounds[key][1] for key in param_bounds],
                           param_opt.n_iteration)
    
    return (xopt,yopt)

def Gillespie_RMMSE(r_diffAB = 0.0328, r_diffBC = 0.0167, alpha_diffAB = 21.10,alpha_diffBC = 46.69, plot = True, plot_RMSE = True):
    
    param_opt.current_iter += 1
    
    print()
    print('Optimization Iteration %s/%s' % (param_opt.current_iter,param_opt.n_iteration))
    
    if plot:
        print("\n=================================================================================") 
        print("========================== Rejection Gillepsie algorithm ========================") 
        print("=================================================================================\n\n")                 
       

    #initialise reactants
    N_init = dict()
    N_init['ESC'] = 100
    N_init['EPI'] = 0
    N_init['NPC'] = 0
    
    #initialise reaction channels
    EE_differentiation = gil.Reaction_channel(param,rate=r_diffAB, shape_param=alpha_diffAB, distribution = param.distribution, name='Differentiation: ESC -> EPI')
    EE_differentiation.reactants = ['ESC']
    EE_differentiation.products = ['EPI']
    
    EN_differentiation = gil.Reaction_channel(param,rate=r_diffBC, shape_param=alpha_diffBC, distribution = param.distribution, name='Differentiation: EPI -> NPC')
    EN_differentiation.reactants = ['EPI']
    EN_differentiation.products = ['NPC']
    reaction_channel_list = [EE_differentiation,EN_differentiation]
    
    #initialise the Gillespie simulation
    G_simul = gil.Gillespie_simulation(N_init,param)
    G_simul.reaction_channel_list = reaction_channel_list
    
    #run multiple Gillespie simulationand average them
    try:
        population = G_simul.run_simulations(param.Tend,verbose = plot)
        if plot:
            G_simul.plot_inter_event_time_distribution(color_list = ['orange','royalblue'])
            G_simul.plot_populations(color_list = ['orange','royalblue', 'gray'])
        
        #load setm cell data
        measured_population = pd.read_csv('Data/Stem_cell_differentiation/stem_cell_data.csv')
        rmse = plot_sim_vs_exp(population, measured_population, reactant_list = ['ESC','EPI','NPC'], cell_line='R1', plot = plot_RMSE)
    
    except:
        rmse = 10
        
    print("RMSE = %.4f \n for r_diffAB = %.4f, r_diffBC = %.4f, \n alpha_diffAB = %.2f, alpha_diffBC = %.2f" % (rmse, r_diffAB, r_diffBC, alpha_diffAB, alpha_diffBC))
    
    
    return rmse
    
    
def plot_results(population,reactant_list, log_scale=False):
    
    """ploting the population"""
    N_simulations = population.shape[0]
    N_reactants = population.shape[2]
    timepoints = population.shape[1]
    time_points = np.linspace(0, param.Tend, timepoints)
    lwm = 3
    plt.figure(figsize = (7,4))
    for ri in range(N_reactants):
        for i in range(N_simulations):
            plt.plot(time_points, population[i,:,ri], 'k-', lw=0.3, alpha=0.2,color=sns.color_palette()[0])
        plt.plot(time_points, population[:,:,ri].mean(axis=0), 'r-', lw=lwm, color=sns.color_palette()[ri+1], label=reactant_list[ri])
    plt.xlabel('Time [min]')
    plt.ylabel('Population')
    plt.legend()
    if log_scale: plt.yscale('log')
    plt.show()
    
    
def plot_sim_vs_exp(population, measured_pop,reactant_list, cell_line = 'R1', xlim=175, filepath=None, plot = True):
    """
    Args:
        cell_line (str): [E14, R1] specify which cell line in the experimental data to use
    """
    
    if plot:
        fig, ax = plt.subplots(1, 3)
        fig.set_size_inches(10,4)
    
        c = [sns.color_palette()[3], sns.color_palette()[2], sns.color_palette()[4]]
        
        timepoints = population.shape[1]
        time_points = np.linspace(0, param.Tend, timepoints)
        
        for ri, reactant in enumerate(reactant_list):
            if reactant != 'Total':
                sns.lineplot(x=time_points, y=population[:,:,ri].mean(axis=0)/np.max(population),
                                 color=c[ri], ax = ax[ri])
                
                # plot markers with CIs
                sns.lineplot(x="time", y="value",  marker='o', linestyle='',
                             color=sns.color_palette()[0], ms = 8,fillstyle='none', 
                             mec=sns.color_palette()[0], mew=1.2,
                             err_style='bars', ci=95, n_boot=1000,
                             data=measured_pop[(measured_pop['L1']==cell_line) & (measured_pop['state']==reactant)], ax = ax[ri])
                ax[ri].set_xlim(0,xlim)
                ax[ri].set(xlabel='Time [h]', ylabel='Probability', title=reactant)
        fig.suptitle('{} cell line'.format(cell_line), fontsize=15)
        fig.tight_layout()
        plt.show()
    
    
    simulated_pop = population.mean(axis=0)/np.max(population)
    rmse = calculate_rmse(simulated_pop, measured_pop, verbose = plot)
    
    return rmse
    

def calculate_rmse(simulated_pop, measured_pop, avg=True, verbose=True, weights=None):
    """Input dataframes need to be in wide format
    
    Args:
        weights (list): [esc, epi, npc] float values on how to weigh rmse average
    
    """
    
    
    
    timepoints = simulated_pop.shape[0]
    time_points = np.linspace(0, param.Tend, timepoints)
    
    simulated_pop_df = pd.DataFrame()
    simulated_pop_df['time'] = time_points
    simulated_pop_df['ESC'] = simulated_pop[:,0]
    simulated_pop_df['EPI'] = simulated_pop[:,1]
    simulated_pop_df['NPC'] = simulated_pop[:,2]
    
    measured_pop_dict = dict()
    for index, row in measured_pop.iterrows():
        if row['L1']=='R1':
            if row['time'] not in measured_pop_dict:
                measured_pop_dict[row['time']] = dict()
            if row['state'] not in measured_pop_dict[row['time']]:
                 measured_pop_dict[row['time']][row['state']] = []
            measured_pop_dict[row['time']][row['state']].append(row['value'])
        
        
    #print(measured_pop_dict)
    time_array = []
    pop_ESC = []
    pop_EPI = []
    pop_NPC = []
    for time in measured_pop_dict:
        time_array.append(time)
        pop_ESC.append(np.mean(measured_pop_dict[time]['ESC']))
        pop_EPI.append(np.mean(measured_pop_dict[time]['EPI']))
        pop_NPC.append(np.mean(measured_pop_dict[time]['NPC']))
    

    measured_pop_df = pd.DataFrame()
    measured_pop_df['time'] = np.array(time_array).astype('float')
    measured_pop_df['ESC'] = np.array(pop_ESC).astype('float')
    measured_pop_df['EPI'] = np.array(pop_EPI).astype('float')
    measured_pop_df['NPC'] = np.array(pop_NPC).astype('float')
    measured_pop_df.index = np.arange(len(measured_pop_df))
    
    # merge to match shape
    merged = pd.merge_asof(measured_pop_df, simulated_pop_df, on='time', direction='nearest', tolerance=0.5, suffixes=('_sim', '_exp')).dropna()
    
    y = merged.loc[:, ['ESC_exp','EPI_exp', 'NPC_exp']]
    yhat = merged.loc[:,['ESC_sim','EPI_sim', 'NPC_sim']]

    if sigmoid_adjusted:
        y_sig = y.transform(lambda x: inv_sigmoid(x))
        yhat_sig = yhat.transform(lambda x: inv_sigmoid(x))
        y = y_sig
        yhat = yhat_sig     

    
    rmse = mean_squared_error(y, yhat, squared=False, multioutput='raw_values')
    if weights:
        rmse_avg = mean_squared_error(y, yhat, squared=False, multioutput=weights)
    else:
        rmse_avg = mean_squared_error(y, yhat, squared=False, multioutput='uniform_average')
    if verbose: print('RMSE:\n ESC = {}\n EPI = {}\n NPC = {}\n Average = {}'.format(rmse[0], rmse[1], rmse[2], rmse_avg)) if verbose else None
    if avg:
        return rmse_avg
    else:
        return rmse
    
    
def sigmoid(y):
  new_y = 1 / (1 + np.exp(-y)) 
  return new_y

def inv_sigmoid(y):
    new_y = np.log(y/(1-y))
    new_y[y > 0.99] = 4.6
    new_y[y < 0.01] = -4.6
    
    return new_y



    
if __name__ == "__main__":
    plt.rcParams.update({'font.size': 14})
    main_optimization()
        
        
