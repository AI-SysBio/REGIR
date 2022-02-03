import numpy as np
import random as rd
import matplotlib.pyplot as plt
import random
import scipy.stats as stats

import REGIR as gil


""" --------------------------------------------------------------------------------------------
3 reactions:
    Goff -> Gon,        Gene activation
    Gon -> Goff,        Gene deactivation
    Gon -> Gon + M,     mRNA production
    M -> 0,             mRNA degradation
"""


class param:
    unit = 'min'
    rd.seed(101)                #length of the simulation, in hours (41 days)
    N_simulations = 200          #The simulation results should be averaged over many trials
    timepoints = 200            #Number of timepoints to record (make surethat this number isnt too big)  
    
    
def main():
    
    r_on = 0.14046623
    r_off = 0.5503625  
    r_prod = 3.53245104 
    r_deg = 0.14551901
    
    rmin = np.min([r_on,r_off,r_prod,r_deg])
    Tend = 10/rmin
    
    param_list = dict()
    param_list['r_deg'] = r_deg
    param_list['r_on'] = r_on
    param_list['r_off'] = r_off
    param_list['r_prod'] = r_prod
    param_list['Tend'] = Tend
    Gillespie_kinetics(param_list, plot = True)   



def Gillespie_kinetics(param_list, plot = True, plot_distrib = True, N_simulations = 200):
    
    if plot:
        print("\n=================================================================================") 
        print("========================== Rejection Gillepsie algorithm ========================") 
        print("=================================================================================\n\n")                 

    #initialise reactants
    Tend = param_list['Tend']
    N_init = dict()
    N_init['Goff'] = 1
    N_init['Gon'] = 0
    N_init['mRNA'] = 0
    
    
    r_on = param_list['r_on']
    r_off = param_list['r_off']
    r_prod = param_list['r_prod'] 
    r_deg = param_list['r_deg'] 
   
    if 'alpha_on' in param_list:
        alpha_on = param_list['alpha_on']
    else:
        alpha_on = 1
    if 'alpha_off' in param_list:
        alpha_off = param_list['alpha_off']
    else:
        alpha_off = 1
    if 'alpha_prod' in param_list:
        alpha_prod = param_list['alpha_prod']
    else:
        alpha_prod = 1
    
    param.N_simulations = N_simulations
    param.Tend = Tend
    #initialise reaction channels
    reaction_channel_list = []
    channel = gil.Reaction_channel(param,rate=r_on, shape_param=alpha_on, distribution = 'Gamma', name = 'Gene activation: Goff -> Gon')
    channel.reactants = ['Goff']
    channel.products = ['Gon']
    reaction_channel_list.append(channel)
    
    channel = gil.Reaction_channel(param,rate=r_off, shape_param=alpha_off, distribution = 'Gamma', name = 'Gene deactivation: Gon -> Goff')
    channel.reactants = ['Gon']
    channel.products = ['Goff']
    reaction_channel_list.append(channel)
    
    channel = gil.Reaction_channel(param,rate=r_prod, shape_param=alpha_prod, distribution = 'Gamma',  name='RNA production: Gon -> Gon + mRNA', transfer_identity = True)
    channel.reactants = ['Gon']
    channel.products = ['Gon','mRNA']
    reaction_channel_list.append(channel)
    
    channel = gil.Reaction_channel(param,rate=r_deg, distribution = 'Exponential',  name='RNA degradation: mRNA -> 0')
    channel.reactants = ['mRNA']
    channel.products = []
    reaction_channel_list.append(channel)
    
    
    
    #initialise the Gillespie simulation
    G_simul = gil.Gillespie_simulation(N_init,param, print_warnings = True)
    G_simul.reaction_channel_list = reaction_channel_list
    if plot:
        print(G_simul)
    
    #run multiple Gillespie simulationand average them
    population = G_simul.run_simulations(Tend, verbose = plot)
    population_per_gene = population/G_simul.reactant_population_init['Goff']
    if plot:
        G_simul.plot_inter_event_time_distribution()
        G_simul.plot_populations(reactant_list = ['Gon','Goff'])
        G_simul.plot_populations(reactant_list = ['mRNA'])
        #plot_distributions(population_per_gene[:,:,2], G_simul.Tend)
    
    mean, CV, Entropy, Fano = get_mRNA_distr_properties(population_per_gene[:,:,2], plot = plot_distrib)
    if plot_distrib: 
        print('   mean %.2f' % mean)
        print('   CV %.2f' % CV)
        print('   Fano %.2f' % Fano)
        print('   Entropy %.2f' % Entropy)
        print()
        
    timepoints = population_per_gene.shape[1]
    mRNApop_steady = population_per_gene[:,int(timepoints/3):,2].reshape(-1)
    
    return mean, CV, Entropy, Fano, mRNApop_steady

def get_mRNA_distr_properties(population, plot = True, pop_is_steady = False):
    #define the second half as the steady state
    if pop_is_steady == False:
        timepoints = population.shape[1]
        pop_steady = population[:,int(timepoints/3):] #steady state is defined as the last third of reaction
        #pop_steady_sorted = np.sort(pop_steady, axis = 1)
        #mean_pop_steady = np.mean(pop_steady_sorted, axis = 0)
    else:
        pop_steady = population
    
    mean_pop_steady = pop_steady.reshape(-1)
    
    plt.figure(figsize = (5,2.4))
    max_count = np.max(mean_pop_steady)
    if max_count < 150:
        bins = np.arange(np.max(mean_pop_steady))
    else:
        bins = 20
    hist = plt.hist(mean_pop_steady, bins = bins, edgecolor='black', color = 'gray', alpha = 0.5, label = 'SSA', density = True)
    if plot: 
        plt.xlabel('mRNA counts')
        plt.ylabel('PDF')
        plt.show()
    else:
        plt.close()
        
        
    std = np.std(mean_pop_steady)
    mean = np.mean(mean_pop_steady)
    if mean > 0:
        CV = std/mean
    else:
        CV = 0
        
    Fano = CV*std
    

    hist = plt.hist(mean_pop_steady, bins = 20, edgecolor='black', color = 'gray', alpha = 0.5, label = 'SSA', density = True)
    plt.clf()
    data = hist[0]
    data = data[data>0]
    data = data/data.sum()
    
    
    LDDP = -(data*np.log(data)).sum()
    #limiting density of discrete points (LDDP) 
    #-> Adjustment to the formula of Claude Shannon for differential entropy.
    #(Differential entropy can be negative -> https://en.wikipedia.org/wiki/Differential_entropy)
    #see https://en.wikipedia.org/wiki/Limiting_density_of_discrete_points
    
    #Deriving Shanon entropy directly from histogram is biased
    #https://stats.stackexchange.com/questions/156235/biased-bootstrap-is-it-okay-to-center-the-ci-around-the-observed-statistic/158683#158683
    #https://stats.stackexchange.com/questions/511516/how-to-derive-the-bias-of-an-entropy-estimate
    
    return mean, CV, LDDP, Fano
    
def plot_distributions(population, Tend):
    """ploting histogram across Gillespie simulations"""
    timepoints = population.shape[1]
    time_points = np.linspace(0, Tend, timepoints)
    Nstep = 1
    for ti in range(Nstep):
        t_index = int((ti+1)*timepoints/Nstep) - 1
        plt.figure(figsize = (5,2.4))
        plt.title('Time %.0f [1/r_deg]' % time_points[t_index])

        pop_plot = population[:,t_index]
        plt.hist(pop_plot, edgecolor='black', color = 'gray', alpha = 0.5, label = 'SSA', density = True)
        plt.xlabel('mRNA counts per gene')
        plt.ylabel('PDF')
    
        var = np.var(pop_plot)
        mean = np.mean(pop_plot)
                
        p = mean / var
        if p < 0:
            p=0.0001
        if p > 1:
            p = 0.9999
        r = p * mean / (1-p)
        print('mRNA', p, r)
                    
        x = np.arange(int(np.max(pop_plot)*1.1))
        pmf = stats.nbinom.pmf(x, r, p)
        plt.plot(x,pmf,lw = 2, color = 'green', label = 'fit')
        plt.legend(prop={'size': 12})
                  
        plt.show()
    

    
    
if __name__ == "__main__":
    plt.rcParams.update({'font.size': 16})
    main()
        
        
