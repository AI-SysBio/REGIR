import numpy as np
import pandas as pd
import random as rd
import matplotlib.pyplot as plt
import seaborn as sns
import random

import REGIR as gil


""" --------------------------------------------------------------------------------------------
2 reactions:
    ESC -> EPI,  differentiation
    EPI -> NPC,  differentiation
    
"""

class param:
    Tend = 170
    unit = 'h'
    rd.seed(101)                #length of the simulation, in hours (41 days)
    N_simulations = 20          #The simulation results should be averaged over many trials
    timepoints = 100            #Number of timepoints to record (make surethat this number isnt too big)
    

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
    
    #run multiple Gillespie simulationand average them
    population = G_simul.run_simulations(param.Tend)
    G_simul.plot_inter_event_time_distribution()
    G_simul.plot_populations()
    
    #load setm cell data
    measured_population = pd.read_csv('Data/Stem_cell_differentiation/stem_cell_data.csv')
    plot_sim_vs_exp(population, measured_population, reactant_list = ['ESC','EPI','NPC'], cell_line='R1')
    
    
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
    
    
def plot_sim_vs_exp(population, population_measured,reactant_list, cell_line = 'R1', xlim=175, filepath=None):
    """
    Args:
        cell_line (str): [E14, R1] specify which cell line in the experimental data to use
    """
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
                         data=population_measured[(population_measured['L1']==cell_line) & (population_measured['state']==reactant)], ax = ax[ri])
            ax[ri].set_xlim(0,xlim)
            ax[ri].set(xlabel='Time [h]', ylabel='Probability', title=reactant)
    fig.suptitle('{} cell line'.format(cell_line), fontsize=15)
    fig.tight_layout()
    



    
if __name__ == "__main__":
    plt.rcParams.update({'font.size': 14})
    main()
        
        
