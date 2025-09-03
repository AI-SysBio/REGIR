import numpy as np
import pandas as pd
import random as rd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from math import log, gamma
from scipy.stats import expon, weibull_min
import os,sys

import lib_REGIR_with_delay as gil


""" --------------------------------------------------------------------------------------------
3 reactions:
    A -> B,  differentiation (exponential)
    B -> C,  differentiation (Weibull, alpha < 1)
    C -> D,  differentiation (Gamma)
"""
class param:
    Tend = 10
    unit = 'h'
    rd.seed(101)                #length of the simulation, in hours (41 days)
    N_simulations = 10         #The simulation results should be averaged over many trials
    timepoints = 100            #Number of timepoints to record (make surethat this number isnt too big)

def main():
    
    print("\n=================================================================================") 
    print("========================= non-Markovian Gillepsie algorithm =======================") 
    print("=================================================================================\n\n")                 
    
    r0 = 1 #same rate for every reactions

    #initialise reactants
    N_init = dict()
    N_init['A'] = 500
    N_init['B'] = 0
    N_init['C'] = 0
    
    #initialise reaction channels
    reaction_channel_list = []
    channel = gil.Reaction_channel(param,rate=r0, distribution = 'Exponential',  name='Exponential Differentiation: A -> B', precompute_delay = False)
    channel.reactants = ['A']
    channel.products = ['B']
    reaction_channel_list.append(channel)
    #No shape parameter 
    
    #channel = gil.Reaction_channel(param,rate=r0/2, shape_param=0.5, distribution = 'Weibull', precompute_delay = False,  name='Gamma Differentiation: B -> C')
    #channel.reactants = ['B']
    #channel.products = ['C']
    #reaction_channel_list.append(channel)
    
    channel = gil.Reaction_channel(param,rate=r0/2, shape_param=0.5, distribution = 'Weibull', precompute_delay = True,  name='Weibull Differentiation: B -> C')
    channel.reactants = ['B']
    channel.products = ['C']
    reaction_channel_list.append(channel)
    
    channel = gil.Reaction_channel(param,rate=r0/2, shape_param=15, distribution = 'Gamma', precompute_delay = True, name='Gamma Differentiation: C -> A')
    channel.reactants = ['C']
    channel.products = ['A']
    reaction_channel_list.append(channel)
    
    
    
    
    
    #initialise the Gillespie simulation
    G_simul = gil.Gillespie_simulation(N_init,param, reaction_channel_list, print_warnings = True)
    
    #run multiple Gillespie simulationand average them
    G_simul.run_simulations(param.Tend, delay_method = 'Anderson')
    #G_simul.run_simulations(param.Tend)
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
        
        
