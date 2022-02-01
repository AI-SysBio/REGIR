import numpy as np
import pandas as pd
import random as rd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from math import log, gamma
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error
from scipy import stats

import REGIR as gil


""" --------------------------------------------------------------------------------------------
3 reactions:
    S -> S+Q, asymmetric division
    S -> Q,   differentiation
    Q -> 0,   apoptosis
    
"""

class param:
    Tend = 400
    unit = 'h'
    rd.seed(100)                #length of the simulation, in hours (41 days)
    N_simulations = 5          #The simulation results should be averaged over many trials
    timepoints = 100            #Number of timepoints to record (make surethat this number isnt too big)
    
    

def main():
    
    print("\n=================================================================================") 
    print("=========================== Rejection Gillepsie algorithm =========================") 
    print("=================================================================================\n\n")                 


    r_div = 1/71.6
    r_diff  = r_div/0.4
    distribution_div = 'LogNormal'
    distribution_diff = 'Gamma'
    alpha_div = 0.1
    
    param.Tend = 4/r_div
    
    alpha_list = np.linspace(1,25,num=13)
    #alpha_list = np.arange(1,20)
    alpha_list = [100] #100 to reproduce paper results
    #distribution_div = 'Gamma'
    #alpha_div = 1
    
    #Exponential
    #distribution_div = 'Exponential'
    #distribution_diff = 'Exponential'
    
    
    CC_data_file = "Data/Bacteria_single_cell_data/CC_data/tableOfTauandKappaInv31C.csv"
    CC_data = pd.read_csv(CC_data_file, names = ['Tau[min]','Kappa inv [min]'])
    wait_times = CC_data['Tau[min]'].to_numpy()
    samples, pdf, x = generate_inter_event_time(rate = r_div, shape_param = alpha_div, distribution = 'lognormal')
    bins = np.linspace(0,np.max(wait_times*1.05), num = 47)
    plt.figure(figsize = (6.2, 4.2))
    plt.hist(wait_times, bins = bins, edgecolor='black', color = 'orange', alpha = 1, density = True, label = 'Single cell data')
    #plt.hist(samples, bins = bins, edgecolor='black', color = sns.color_palette()[4], alpha = 0.5, density = True, label = 'Simulated')
    plt.plot(x, pdf, color = 'black', label = 'Fitted (lognormal)')
    plt.legend(loc = 'upper left', prop = {'size':14})
    plt.xlabel('Time to division [min]')
    plt.ylabel('PDF')
    plt.xlim(0,110)
    plt.show()    
    

    
    mean_scores = []
    std_scores = []
    
    for i, alphai in enumerate(alpha_list):
        
        alpha_diff = alphai

        #initialise reactants
        N_init = dict()
        N_init['S'] = 100
        N_init['Q'] = 0
        
        #initialise reaction channels
        S_division = gil.Reaction_channel(param,rate=r_div, shape_param=alpha_div, distribution = distribution_div, name='S -> S+Q')
        S_division.reactants = ['S']
        S_division.products = ['S', 'Q']
        
        SQ_differentiation = gil.Reaction_channel(param,rate=r_diff, shape_param=alpha_diff, distribution = distribution_diff, name='Q -> S')
        SQ_differentiation.reactants = ['Q']
        SQ_differentiation.products = ['S']
    
        reaction_channel_list = [S_division,SQ_differentiation]
        
        #initialise the Gillespie simulation
        G_simul = gil.Gillespie_simulation(N_init,param)
        G_simul.reaction_channel_list = reaction_channel_list
        
        #run multiple Gillespie simulationand average them
        population = G_simul.run_simulations(param.Tend)
        G_simul.plot_inter_event_time_distribution(color_list = ['orange','royalblue'])
        G_simul.plot_populations(color_list = ['orange','royalblue', 'gray'], figsize = (5,4))
        
        mean_score, std_score = quantify_oscillations(population[:,:,1])
        mean_scores.append(mean_score)
        std_scores.append(std_score)
        
    mean_scores = np.array(mean_scores)
    std_scores = np.array(std_scores)
    
    if len(alpha_list) > 2:
        plt.plot(alpha_list,mean_scores,'o-', color = 'black', lw = 2)
        plt.fill_between(alpha_list,mean_scores - std_scores, mean_scores + std_scores, color = 'black', alpha = 0.3)
        plt.xlabel('Cell differentiation shape parameter')
        plt.ylabel('Oscillation Amplitude')
        plt.xticks([1,5,10,15,20,25],[1,5,10,15,20,25])
        plt.show()
        
        
    
    
def plot_results(population,reactant_list, log_scale=False):
    
    """ploting the population"""
    N_simulations = population.shape[0]
    N_reactants = population.shape[2]
    timepoints = population.shape[1]
    time_points = np.linspace(0, param.Tend, timepoints)
    lwm = 3
    plt.figure(figsize = (5,4))
    plt.rcParams.update({'font.size': 16})
    color_list = ['orange','royalblue', 'gray']
    for ri in range(N_reactants):
        for i in range(N_simulations):
            plt.plot(time_points, population[i,:,ri], 'k-', lw=0.3, alpha=0.1,color=sns.color_palette()[0])
        plt.plot(time_points, population[:,:,ri].mean(axis=0), 'r-', lw=lwm, color=color_list[ri], label=reactant_list[ri])
    plt.xlabel('Time [min]')
    plt.ylabel('Population')
    plt.legend()
    if log_scale: plt.yscale('log')
    plt.show()
    
def quantify_oscillations(population_Q, Verbose = False):
    #fit exp function
    def exp_func(x,k,a):
        return k*np.exp(a*x)
    
    timepoints = population_Q.shape[1]
    Nsim = population_Q.shape[0]
    xdata = np.linspace(0, param.Tend, timepoints)
    Scores = []
    for si in range(Nsim):
        ydata = population_Q[si,:]
        popt, pcov = curve_fit(exp_func, xdata, ydata, bounds=(0.000001, 100), p0=[2,8/100])
        yfit = exp_func(xdata, *popt)
        
        if Verbose:
            plt.plot(xdata, ydata, 'b-', label='data')
            plt.plot(xdata, yfit, 'g--')
            plt.xlabel('Time')
            plt.ylabel('Population')
            plt.show()
    
        Score = mean_squared_error(ydata, yfit)
        Scores.append(Score)
    
    return np.mean(Scores), np.std(Scores)



def generate_inter_event_time(rate=1, shape_param=6, distribution = 'Gamma', n = 10000):
    
    if distribution.lower() in ['exponential', 'exp']:
        scale = 1/rate
        samples = stats.expon.rvs(loc=0, scale=scale, size=n)    
        pdf = stats.expon.pdf(x=np.linspace(0,110,num = 1000), loc=0, scale=scale)    
    
    elif distribution.lower() in ['gamma','gam']:
        alpha = shape_param
        beta = alpha*rate
        scale = 1 / beta
        samples = stats.gamma.rvs(alpha, loc=0, scale=scale, size=n) 
        pdf = stats.gamma.pdf(x=np.linspace(0,110,num = 1000), a=alpha, loc=0, scale=scale)
        
    elif distribution.lower() in ['weibull','weib']:
        k = shape_param
        r0 = rate
        beta = k * (r0 * gamma((k + 1)/k))**(k)
        lam = np.power(k/beta, 1/k) # scale
        samples = stats.weibull_min.rvs(k, loc=0, scale=lam, size=n)    
        pdf = stats.weibull_min.pdf(np.linspace(0,110,num = 1000), k, loc=0, scale=scale)
        
    elif distribution.lower() in ['gaussian', 'normal', 'norm']:
        samples = stats.norm.rvs(loc=1/rate, scale=shape_param, size=n)
        pdf = stats.norm.pdf(x=np.linspace(0,110,num = 1000), loc=1/rate, scale=shape_param)
        
    elif distribution.lower() in ['lognormal', 'lognorm']:
        sigma0 = 1/rate*shape_param
        mu0 = 1/rate
        mu = log(mu0**2 / np.sqrt(mu0**2 + sigma0**2))
        sigma = np.sqrt(log(1+ sigma0**2/mu0**2))
        samples = stats.lognorm.rvs(s= sigma, loc=0, scale=np.exp(mu), size=n)
        pdf =  stats.lognorm.pdf(x=np.linspace(0,110,num = 1000), s= sigma, loc=0, scale=np.exp(mu))
        
        
        """
        sigma0 = 1/rate*alpha
        mu0 = 1/rate
        mu = log(mu0**2 / sqrt(mu0**2 + sigma0**2))
        sigma = sqrt(log(1+ sigma0**2/mu0**2))
        if t == 0:
            rate = 0
        else:
            pdf = 1/(t*sigma*sqrt(2*pi)) * exp(-(1/2) * (log(t)-mu)**2/(sigma**2))
            cdf = 1/2*(1+erf((log(t)-mu)/(sigma*sqrt(2))))
            rate = pdf/(1-cdf)
        """
        
        
    return samples, pdf, np.linspace(0,110,num = 1000)
    
    

    
if __name__ == "__main__":
    plt.rcParams.update({'font.size': 16})
    main()
        
        
