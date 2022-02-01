import numpy as np
import pandas as pd
import random as rd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from math import log, gamma
from scipy.stats import expon, weibull_min
import os,sys
from datetime import datetime

import lib_nMGA as gil_nMGA
import lib_REGIR as gil_REGIR


""" --------------------------------------------------------------------------------------------
1 reactions:
    A -> B,  differentiation (Gamma)
    
"""
class param:
    Tend = 5
    unit = 'h'
    rd.seed(101)                #length of the simulation, in hours (41 days)
    N_simulations = 10         #The simulation results should be averaged over many trials
    timepoints = 100            #Number of timepoints to record (make surethat this number isnt too big)

def main():
    N_list = np.logspace(3,18,num = 16, base=2).astype(int)
    
    recompute_time_complexity = False
    if recompute_time_complexity:
        compute_time_complexity(N_list)
    
    time_Exp = np.load('time_EXP.npy')
    time_REGIR = np.load('time_REGIR.npy')
    time_nMGA = np.load('time_nMGA.npy')
    
    #time_nMGA = np.array([0.0800621000,0.367517000,1.28012000,4.93146800,17.7352620,97.6663530, 3.53343002e+02, 1.30758699e+03, 5.16693849e+03, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
    #time_REGIR = np.array([1.38576000e-01, 2.79990000e-01, 5.51963000e-01, 1.188415000e+00, 2.27323500e+00, 5.2836460e+00, 1.09630070e+01, 1.67750130e+01,3.28149910e+01, 7.21816420e+01, 1.26181540e+02, 2.57781622e+02, 4.71099330e+02, 1.21256787e+03, 2.20974554e+03, 4.47487488e+03])
    #time_Exp = np.array([5.36930000e-02, 1.065360000e-01, 2.09594000e-01, 4.08914000e-01, 7.82378000e-01, 1.67291500e+00, 3.59508100e+00, 6.93443900e+00, 1.26704890e+01, 2.35766060e+01, 4.58852610e+01, 8.79955970e+01, 1.88073875e+02, 3.83659120e+02, 6.59597018e+02, 1.45400972e+03])
    
    time_nMGA = np.append(time_nMGA,np.nan)
    time_nMGA = np.append(time_nMGA,np.nan)
    time_nMGA = np.append(time_nMGA,np.nan)
    time_REGIR = np.append(time_REGIR,np.nan)
    time_REGIR = np.append(time_REGIR,np.nan)
    time_REGIR = np.append(time_REGIR,np.nan)
    time_Exp = np.append(time_Exp,np.nan)
    time_Exp = np.append(time_Exp,np.nan)
    time_Exp = np.append(time_Exp,np.nan)
     
    
    N_list = np.logspace(3,21,num = 19, base=2).astype(int)
    
    
    fitted_Exp = fit_log_line(N_list,time_Exp)
    fitted_REGIR = fit_log_line(N_list,time_REGIR)
    fitted_nMGA = fit_log_line(N_list,time_nMGA)

    plt.figure(figsize = (8,4))
    plt.scatter(N_list, time_Exp/10, s=60, color = 'blue', label = 'Gillespie (exponential)')
    plt.plot(N_list, fitted_Exp/10, lw=2, color = 'blue')
    plt.scatter(N_list, time_REGIR/10, s=60, color = 'green', label = r'REGIR (gamma, $\alpha=6$)')
    plt.plot(N_list, fitted_REGIR/10, lw=2, color = 'green')
    plt.scatter(N_list, time_nMGA/10, s=60, color = 'red', label = r'nMGA (gamma, $\alpha=6$)')
    plt.plot(N_list, fitted_nMGA/10, lw=2, color = 'red')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Number of reactants')
    plt.ylabel('Av. time / simulation [s]')
    plt.ylim(2e-3,2e3)
    plt.xlim(5e0,1e6)
    plt.legend(loc = 'lower right', prop={"size":14})
    plt.show()
    
def fit_log_line(x0,y0):
    
    x1 = x0[~np.isnan(y0)]
    y1 = y0[~np.isnan(y0)]
    
    x = np.log(x1)
    y = np.log(y1)
    m, b = np. polyfit(x, y, 1)
    print("  Fitted slope %.2f" % m)
    
    fitted_y = np.exp(m*np.log(x0) + b)
    return fitted_y
    
    
        

def compute_time_complexity(N_list):
    
    print("\n=================================================================================") 
    print("========================= non-Markovian Gillepsie algorithm =======================") 
    print("=================================================================================\n\n")                 
    
    r0 = 1 #same rate for every reactions
    
    time_REGIR = []
    time_nMGA = []
    time_Exp = []
    
    for ni,N in enumerate(N_list):

        #initialise reactants
        N_init = dict()
        N_init['A'] = N
        
        
        """ Gillespie Exp algorithm benchmarking"""
        #initialise reaction channels
        reaction_channel_list = []
        channel = gil_REGIR.Reaction_channel(param,rate=r0, distribution = 'Exponential')
        channel.reactants = ['A']
        channel.products = ['A']
        reaction_channel_list.append(channel)


        #initialise the Gillespie simulation
        G_simul = gil_REGIR.Gillespie_simulation(N_init,param, print_warnings = True)
        G_simul.reaction_channel_list = reaction_channel_list
        
        #run multiple Gillespie simulationand average them
        t_start = datetime.now()
        G_simul.run_simulations(param.Tend, verbose=False)
        t_end = datetime.now()
        t_elapsed = t_end - t_start
        
        time_Exp.append(t_elapsed.total_seconds())
        
        #G_simul.plot_inter_event_time_distribution()
        #G_simul.plot_populations(figsize = (8,4))
    
        if N < 3000:
            """ nMGA algorithm benchmarking"""
            #initialise reaction channels
            reaction_channel_list = []
            channel = gil_nMGA.Reaction_channel(param,rate=r0, shape_param=6, distribution = 'Gamma')
            channel.reactants = ['A']
            channel.products = ['A']
            reaction_channel_list.append(channel)
    
    
            #initialise the Gillespie simulation
            G_simul = gil_nMGA.Gillespie_simulation(N_init,param, print_warnings = True)
            G_simul.reaction_channel_list = reaction_channel_list
            
            #run multiple Gillespie simulationand average them
            t_start = datetime.now()
            G_simul.run_simulations(param.Tend, verbose=False)
            t_end = datetime.now()
            t_elapsed = t_end - t_start
            
            time_nMGA.append(t_elapsed.total_seconds())
            
            #G_simul.plot_inter_event_time_distribution()
            #G_simul.plot_populations(figsize = (8,4))
        else:
            time_nMGA.append(np.nan)
        
        
        
        """ REGIR algorithm benchmarking"""
            #initialise reaction channels
        reaction_channel_list = []
        channel = gil_REGIR.Reaction_channel(param,rate=r0, shape_param=6, distribution = 'Gamma')
        channel.reactants = ['A']
        channel.products = ['A']
        reaction_channel_list.append(channel)


        #initialise the Gillespie simulation
        G_simul = gil_REGIR.Gillespie_simulation(N_init,param, print_warnings = True)
        G_simul.reaction_channel_list = reaction_channel_list
        
        #run multiple Gillespie simulationand average them
        t_start = datetime.now()
        G_simul.run_simulations(param.Tend, verbose=False)
        t_end = datetime.now()
        t_elapsed = t_end - t_start
        
        time_REGIR.append(t_elapsed.total_seconds())
        
        #G_simul.plot_inter_event_time_distribution()
        #G_simul.plot_populations(figsize = (8,4))
        
        
        
        print("For N=%s, nMGA elapsed time was %s seconds for %s simulations" % (N,time_nMGA[-1],param.N_simulations))
        print("For N=%s, REGIR elapsed time was %s seconds for %s simulations" % (N,time_REGIR[-1],param.N_simulations))
        print("For N=%s, Gillespie elapsed time was %s seconds for %s simulations" % (N,time_Exp[-1],param.N_simulations))
        print()
        
        
        
        
    np.save('time_EXP.npy', time_Exp)
    np.save('time_REGIR.npy', time_REGIR)
    np.save('time_nMGA.npy', time_nMGA)
    

    
    
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
    plt.rcParams.update({'font.size': 16})
    main()
        
        
