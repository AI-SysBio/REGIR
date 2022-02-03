import numpy as np
import random as rd
import matplotlib.pyplot as plt
import seaborn as sns
import random

import REGIR as gil


""" --------------------------------------------------------------------------------------------
6 reactions:
    A -> B,  differentiation (exponential)
    B -> C,  differentiation (Weibull)
    C -> D,  differentiation (Gaussian)
    D -> E,  differentiation (Gamma)
    E -> F,  differentiation (Cauchy)
    F -> A,  differentiation (LogNormal)
    
"""
class param:
    Tend = 10
    unit = 'h'
    rd.seed(101)                #length of the simulation, in hours (41 days)
    N_simulations = 20          #The simulation results should be averaged over many trials
    timepoints = 100            #Number of timepoints to record (make surethat this number isnt too big)

def main():
    
    print("\n=================================================================================") 
    print("========================= non-Markovian Gillepsie algorithm =======================") 
    print("=================================================================================\n\n")                 
    
    r0 = 1 #same rate for every reactions

    #initialise reactants
    N_init = dict()
    N_init['A'] = 100
    N_init['B'] = 0
    N_init['C'] = 0
    N_init['D'] = 0
    N_init['E'] = 0
    N_init['F'] = 0
    
    #initialise reaction channels
    reaction_channel_list = []
    channel = gil.Reaction_channel(param,rate=r0, distribution = 'Exponential',  name='Exponential Differentiation: A -> B')
    channel.reactants = ['A']
    channel.products = ['B']
    reaction_channel_list.append(channel)
    #No shape parameter 
    
    channel = gil.Reaction_channel(param,rate=3*r0, shape_param=2.5, distribution = 'Weibull',  name='Weibull Differentiation: B -> C')
    channel.reactants = ['B']
    channel.products = ['C']
    reaction_channel_list.append(channel)
    #shape parameter rate = t^(alpha-1)
    
    channel = gil.Reaction_channel(param,rate=r0/2, shape_param=15, distribution = 'Gamma',  name='Gamma Differentiation: C -> D')
    channel.reactants = ['C']
    channel.products = ['D']
    reaction_channel_list.append(channel)
    #shape parameter rate = k in gamma parametrization
    
    channel = gil.Reaction_channel(param,rate=10*r0, shape_param=0.2, distribution = 'Gaussian',  name='Gaussian Differentiation: D -> E')
    channel.reactants = ['D']
    channel.products = ['E']
    reaction_channel_list.append(channel)
    #shape parameter = standard deviation sigma
    
    channel = gil.Reaction_channel(param,rate=r0, shape_param=0.2, distribution = 'Cauchy',  name='Cauchy Differentiation: E -> F')
    channel.reactants = ['E']
    channel.products = ['F']
    reaction_channel_list.append(channel)
    #shape parameter = gamma
    
    
    channel = gil.Reaction_channel(param,rate=3*r0, shape_param=5, distribution = 'Weibull',  name='LogNormal Differentiation: F -> A')
    channel.reactants = ['F']
    channel.products = ['A']
    reaction_channel_list.append(channel)
    #shape parameter = sigma
    
    
    
    
    #initialise the Gillespie simulation
    G_simul = gil.Gillespie_simulation(N_init,param, print_warnings = True)
    G_simul.reaction_channel_list = reaction_channel_list
    print(G_simul)
    
    #run multiple Gillespie simulationand average them
    G_simul.run_simulations(param.Tend)
    G_simul.plot_inter_event_time_distribution()
    G_simul.plot_populations(figsize = (8,4))
    
    
def plot_results(population,reactant_list, log_scale=False):
    
    """ploting the population"""
    N_simulations = population.shape[0]
    N_reactants = population.shape[2]
    timepoints = population.shape[1]
    time_points = np.linspace(0, param.Tend, timepoints)
    lwm = 3
    plt.figure(figsize = (8,4))
    for ri in range(N_reactants):
        for i in range(N_simulations):
            plt.plot(time_points, population[i,:,ri], 'k-', lw=0.3, alpha=0.2,color=sns.color_palette()[0])
        plt.plot(time_points, population[:,:,ri].mean(axis=0), 'r-', lw=lwm, color=sns.color_palette()[ri+1], label=reactant_list[ri])
    plt.xlabel('Time [h]')
    plt.ylabel('Population')
    plt.legend(prop={'size': 12})
    if log_scale: plt.yscale('log')
    plt.show()
    
    
if __name__ == "__main__":
    plt.rcParams.update({'font.size': 14})
    main()
        
        
