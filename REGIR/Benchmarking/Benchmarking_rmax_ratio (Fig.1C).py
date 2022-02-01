import numpy as np
import pandas as pd
import random as rd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from math import log, gamma
from scipy import stats
from scipy.optimize import curve_fit
import os,sys

import lib_REGIR as gil_REGIR


""" --------------------------------------------------------------------------------------------
1 reactions:
    A -> B,  differentiation (Gamma)
    
"""
class param:
    Tend = 5
    unit = 'h'
    rd.seed(101)                #length of the simulation, in hours (41 days)
    N_simulations = 10        #The simulation results should be averaged over many trials
    timepoints = 100            #Number of timepoints to record (make surethat this number isnt too big)


def main():
    
    if not os.path.isdir("Computed"):
        os.makedirs("Computed")    

    shape_param_list = np.linspace(1,15,num = 10)
    shape_param_list_2 = np.logspace(-1,0,base=10,num=10)[::-1]
    
    recompute_time_complexity = False
    if recompute_time_complexity:
        compute_rmax_ratio(shape_param_list,shape_param_list_2)
    
    ratio_weibull = np.load('Computed/Ratio_weibull.npy')   
    ratio_gamma = np.load('Computed/Ratio_gamma.npy')
    ratio_normal = np.load('Computed/Ratio_normal.npy')
    ratio_lognormal = np.load('Computed/Ratio_lognormal.npy')
    
    
    fitted_weibull = fit_line(shape_param_list,ratio_weibull)
    #fitted_gamma = fit_line(shape_param_list,ratio_gamma)
    

    plt.figure(figsize = (9,4))
    plt.axhline(y=0, color = 'black', linestyle = '--', label = 'Exponential', alpha = 0.5)
    plt.plot(shape_param_list, ratio_weibull, 'o-', markersize=8, color = 'purple', linestyle='dashed', label = r'Weibull ($\alpha$)')
    plt.plot(shape_param_list, fitted_weibull, lw=1, color = 'purple', alpha = 0.4)
    plt.plot(shape_param_list, ratio_gamma, 'o-', markersize=8, linestyle='dashed', color = 'orange', label = r'Gamma ($\alpha$)')
    #plt.plot(shape_param_list, fitted_gamma, lw=2, color = 'orange', linestyle = '--')
    plt.plot(shape_param_list, ratio_normal, 'o-', markersize=8, linestyle='dashed', color = 'green', label = r'Normal ($\sigma$)')
    plt.plot(shape_param_list, ratio_lognormal, 'o-', markersize=8, linestyle='dashed', color = sns.color_palette()[3], label = r'LogNormal ($\sigma$)')
    plt.xlabel(r'Shape parameter')
    plt.ylabel('Rejected/accepted ratio')
    plt.legend(loc = 'upper left', prop={"size":13})
    plt.xticks(ticks = np.arange(8)*2+1, labels = np.arange(8)*2+1)
    plt.xticks(ticks = np.arange(8)*2+1, labels = np.logspace(-1,0,base=10,num=8)[::-1].round(2))
    plt.xlim(0.3,15.7)
    plt.show()
    
def fit_line(x0,y0):
    
    x0 = np.array(x0)
    y0 = np.array(y0)
    
    def fit_func(x, a):
        # Curve fitting function
        return a*(x-1)
    
    m = curve_fit(fit_func, x0, y0)[0][0]
    print(m)
    fitted_y = m*(x0-1)
    return fitted_y


def compute_rmax_ratio(shape_param_list,shape_param_list_2):
    
    print("\n=================================================================================") 
    print("======================== non-Markovian Gillepsie algorithm ======================") 
    print("=================================================================================\n\n")                 
    
    r0 = 1 #same rate for every reactions
    
    ratio_gamma = []
    ratio_weibull = []
    ratio_normal = []
    ratio_lognormal = []
    

    #initialise reactants
    N_init = dict()
    N_init['A'] = 500
    
    for ai, shape_param in enumerate(shape_param_list):
        
        alpha = shape_param_list[ai]
        sigma = shape_param_list_2[ai]

        """ rejected/accepted ratio for Gamma distribution """
        reaction_channel_list = []
        channel = gil_REGIR.Reaction_channel(param,rate=r0, shape_param=alpha, distribution = 'Gamma')
        channel.reactants = ['A']
        channel.products = ['A']
        reaction_channel_list.append(channel)

        #initialise the Gillespie simulation
        G_simul = gil_REGIR.Gillespie_simulation(N_init,param, print_warnings = True, min_ratio = 1)
        G_simul.reaction_channel_list = reaction_channel_list
        G_simul.run_simulations(param.Tend, verbose=False)
        
        print(reaction_channel_list[0].number_of_rejected_reactions)
        print(reaction_channel_list[0].number_of_accepted_reactions)
        
        ratio = reaction_channel_list[0].number_of_rejected_reactions/reaction_channel_list[0].number_of_accepted_reactions
        ratio_gamma.append(ratio)
        G_simul.plot_inter_event_time_distribution()
        
        
        
        """ rejected/accepted ratio for Weibull distribution """
        reaction_channel_list = []
        channel = gil_REGIR.Reaction_channel(param,rate=r0, shape_param=alpha, distribution = 'Weibull')
        channel.reactants = ['A']
        channel.products = ['A']
        reaction_channel_list.append(channel)

        #initialise the Gillespie simulation
        G_simul = gil_REGIR.Gillespie_simulation(N_init,param, print_warnings = True, min_ratio = 0) #no approximation correction
        G_simul.reaction_channel_list = reaction_channel_list
        G_simul.run_simulations(param.Tend, verbose=False)
        
        print(reaction_channel_list[0].number_of_rejected_reactions)
        print(reaction_channel_list[0].number_of_accepted_reactions)        

        ratio = reaction_channel_list[0].number_of_rejected_reactions/reaction_channel_list[0].number_of_accepted_reactions
        ratio_weibull.append(ratio)
        G_simul.plot_inter_event_time_distribution()
        
        
        """ rejected/accepted ratio for Normal distribution """
        reaction_channel_list = []
        channel = gil_REGIR.Reaction_channel(param,rate=r0, shape_param=sigma, distribution = 'normal')
        channel.reactants = ['A']
        channel.products = ['A']
        reaction_channel_list.append(channel)

        #initialise the Gillespie simulation
        G_simul = gil_REGIR.Gillespie_simulation(N_init,param, print_warnings = True, min_ratio = 0) #no approximation correction
        G_simul.reaction_channel_list = reaction_channel_list
        G_simul.run_simulations(param.Tend, verbose=False)
        
        print(reaction_channel_list[0].number_of_rejected_reactions)
        print(reaction_channel_list[0].number_of_accepted_reactions)        

        ratio = reaction_channel_list[0].number_of_rejected_reactions/reaction_channel_list[0].number_of_accepted_reactions
        ratio_normal.append(ratio)
        G_simul.plot_inter_event_time_distribution()
        
        
        
        
        """ rejected/accepted ratio for LogNormal distribution """
        reaction_channel_list = []
        channel = gil_REGIR.Reaction_channel(param,rate=r0, shape_param=sigma, distribution = 'lognormal')
        channel.reactants = ['A']
        channel.products = ['A']
        reaction_channel_list.append(channel)

        #initialise the Gillespie simulation
        G_simul = gil_REGIR.Gillespie_simulation(N_init,param, print_warnings = True, min_ratio = 0) #no approximation correction
        G_simul.reaction_channel_list = reaction_channel_list
        G_simul.run_simulations(param.Tend, verbose=False)
        
        print(reaction_channel_list[0].number_of_rejected_reactions)
        print(reaction_channel_list[0].number_of_accepted_reactions)        

        ratio = reaction_channel_list[0].number_of_rejected_reactions/reaction_channel_list[0].number_of_accepted_reactions
        ratio_lognormal.append(ratio)
        G_simul.plot_inter_event_time_distribution()
        
        

        print("For alpha=%.2f, Gamma ratio is %.2f" % (alpha,ratio_gamma[-1]))
        print("For alpha=%.2f, Weibull ratio is %.2f" % (alpha,ratio_weibull[-1]))
        print("For sigma=%.2f, Normal ratio is %.2f" % (sigma,ratio_normal[-1]))
        print("For sigma=%.2f, LogNormal ratio is %.2f" % (sigma,ratio_lognormal[-1]))
        print()
        

    np.save('Ratio_weibull.npy', ratio_weibull)        
    np.save('Ratio_gamma.npy', ratio_gamma)
    np.save('Ratio_normal.npy', ratio_normal)        
    np.save('Ratio_lognormal.npy', ratio_lognormal)
    

    
    
if __name__ == "__main__":
    plt.rcParams.update({'font.size': 16})
    main()
        
        
