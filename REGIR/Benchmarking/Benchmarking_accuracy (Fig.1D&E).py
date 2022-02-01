import numpy as np
import pandas as pd
import random as rd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from math import log, gamma
from scipy import stats
import os,sys
from datetime import datetime

from pyemd import emd_samples
import REGIR as gil_REGIR


""" --------------------------------------------------------------------------------------------
1 reactions:
    A -> B,  differentiation (Weibull)
    
"""

shape_param = 6

class param:
    Tend = 5
    unit = 'h'
    rd.seed(110)                #length of the simulation, in hours (41 days)
    N_simulations = 10        #The simulation results should be averaged over many trials
    timepoints = 100            #Number of timepoints to record (make surethat this number isnt too big)



def main():
    
    if not os.path.isdir("Computed"):
        os.makedirs("Computed")
    
    """
    #quick check of the abckground intrisic to EMD
    n = 10000
    r0 = 1
    y1 = generate_inter_event_time(rate=r0, shape_param=shape_param, distribution = 'gamma', n=n)
    y2 = generate_inter_event_time(rate=r0, shape_param=shape_param, distribution = 'gamma', n=n)
    EMD,std = emd_av(y1, y2)  
    print(EMD*100)
    """

    N_list = np.logspace(0,7,num = 8, base=2).astype(int)

    recompute_accuracy = False
    if recompute_accuracy:
        compute_EMD(N_list)
    
    EMDrej = np.load('Computed/EMDrej.npy')   
    EMD1 = np.load('Computed/EMD1.npy')

    plt.figure(figsize = (8,3.5))
    plt.plot(N_list, EMD1, 'o-', markersize=8, linestyle='dashed', color = sns.color_palette()[4], label = r'REGIR ($\lambda_{max} \geq \lambda_0$)')
    plt.plot(N_list, EMDrej, 'o-', markersize=8, linestyle='dashed', color = sns.color_palette()[5], label = r'REGIR adjusted ($f=30$)')
    plt.fill_between(np.logspace(-1,8,num = 2, base=2), 0,1.5, color = 'black', alpha = 0.2, label = 'Statistically insignificant')
    plt.xlabel(r'Number of reactants')
    plt.ylabel('REGIR error (EMD) [%]')
    plt.legend(loc = 'upper right', prop={"size":14})
    plt.xlim(8.5e-1,1e2)
    plt.xscale('log')
    plt.show()


def compute_EMD(N_list):
    
    print("\n=================================================================================") 
    print("======================== non-Markovian Gillepsie algorithm ======================") 
    print("=================================================================================\n\n")                 
    
    r0 = 1 #same rate for every reactions
    y_true = generate_inter_event_time(rate=r0, shape_param=shape_param, distribution = 'gamma')
    EMDrej = []
    EMD1 = []
    #EMD2 = []
    
    for ni,N in enumerate(N_list):
        
        param.N_simulations = int(10000/N)+1
        #param.N_simulations = 100


        #initialise reactants
        N_init = dict()
        N_init['A'] = N
        
        
        """ Rejection first order approximation """
        reaction_channel_list = []
        channel = gil_REGIR.Reaction_channel(param,rate=r0, shape_param=shape_param, distribution = 'gamma')
        channel.reactants = ['A']
        channel.products = []
        reaction_channel_list.append(channel)

        #initialise the Gillespie simulation
        G_simul = gil_REGIR.Gillespie_simulation(N_init,param, print_warnings = True, min_ratio = 30)
        G_simul.reaction_channel_list = reaction_channel_list
        G_simul.run_simulations(param.Tend, verbose=False)
        
        y_sim = np.array(G_simul.reaction_channel_list[0].wait_times)
        EMD,std = emd_av(y_sim, y_true)
        normed_EMD = EMD*r0
        EMDrej.append(normed_EMD*100)
        
        #G_simul.plot_inter_event_time_distribution(plot_fitted = False, theory_color = 'black', bins = np.linspace(0,5,40))
        #G_simul.plot_populations(figsize = (8,4))
        
        
        
        """ First order approximation"""
        reaction_channel_list = []
        channel = gil_REGIR.Reaction_channel(param,rate=r0, shape_param=shape_param, distribution = 'gamma')
        channel.reactants = ['A']
        channel.products = []
        reaction_channel_list.append(channel)

        #initialise the Gillespie simulation
        G_simul = gil_REGIR.Gillespie_simulation(N_init,param, print_warnings = True, min_ratio = 0) #no approximation correction
        G_simul.reaction_channel_list = reaction_channel_list
        G_simul.run_simulations(param.Tend, verbose=False)
        
        y_sim = np.array(G_simul.reaction_channel_list[0].wait_times)
        EMD,std = emd_av(y_sim, y_true)
        normed_EMD = EMD*r0
        EMD1.append(normed_EMD*100)
        
        G_simul.plot_inter_event_time_distribution(plot_fitted = False, theory_color = 'black', bins = np.linspace(0,5,40))
        
        
        print("For N=%s, REGIR rej error is %s" % (N,EMDrej[-1]))
        print("For N=%s, REGIR raw error is %s" % (N,EMD1[-1]))
        print()
        

    np.save('EMDrej.npy', EMDrej)        
    np.save('EMD1.npy', EMD1)
    
    
def emd_av(y1,y2):
    
    print(len(y1))
    print(len(y2))
    return emd_samples(y1,y2),0
    
    
def generate_inter_event_time(rate=1, shape_param=6, distribution = 'Gamma', n = 10000):
    
    
    if distribution.lower() in ['exponential', 'exp']:
        scale = 1/rate
        y_true = stats.expon.rvs(loc=0, scale=scale, size=n)       
    
    elif distribution.lower() in ['gamma','gam']:
        alpha = shape_param
        beta = alpha*rate
        scale = 1 / beta
        y_true = stats.gamma.rvs(alpha, loc=0, scale=scale, size=n)  
        
    elif distribution.lower() in ['weibull','weib']:
        k = shape_param
        r0 = rate
        beta = k * (r0 * gamma((k + 1)/k))**(k)
        lam = np.power(k/beta, 1/k) # scale
        y_true = stats.weibull_min.rvs(k, loc=0, scale=lam, size=n)    
        
    elif distribution.lower() in ['gaussian', 'normal', 'norm']:
        y_true = stats.norm.rvs(loc=1/rate, scale=shape_param, size=n) 
        
        
    return y_true
    

    
    
if __name__ == "__main__":
    plt.rcParams.update({'font.size': 14})
    main()
        
        
