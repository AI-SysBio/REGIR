import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stat
from math import log, gamma, exp, factorial, pi, sqrt, erf, atan
from scipy.special import gammainc

from pyemd import emd_samples

def main():
    
    #single cell division experimental data from https://jun.ucsd.edu/repository.php
    experimental_waiting_times = []
    B_subtilis_data_file = 'Data/Bacteria_single_cell_data/B-subtilis.xlsx'
    B_data = pd.read_excel(B_subtilis_data_file, sheet_name=3)
    wait_times = B_data['generation time (minute)'].to_numpy()
    experimental_waiting_times.append(wait_times)
    
    B_C_data_file = 'Data/Bacteria_single_cell_data/E-coli.xlsx'
    BC_data = pd.read_excel(B_C_data_file, sheet_name=3)
    wait_times = BC_data['generation time (minute)'].to_numpy()
    wait_times = wait_times[wait_times<200]
    experimental_waiting_times.append(wait_times)
    
    CC_data_file = "Data/Bacteria_single_cell_data/CC_data/tableOfTauandKappaInv31C.csv"
    CC_data = pd.read_csv(CC_data_file, names = ['Tau[min]','Kappa inv [min]'])
    wait_times = CC_data['Tau[min]'].to_numpy()
    experimental_waiting_times.append(wait_times)

    distribution_list = {'Gaussian':stat.norm, 'Weibull':stat.weibull_min, 'Gamma':stat.gamma, 'LogNormal':stat.lognorm}


    for wi,wait_times in enumerate(experimental_waiting_times):
    
        fig1, ax1 = plt.subplots()
        ax1.set_xlabel('Time [min]')
        ax1.set_ylabel('PDF')
        fig2, ax2 = plt.subplots()
        ax2.set_xlabel('Time [min]')
        ax2.set_ylabel('CDF')
        fig3, ax3 = plt.subplots()
        ax3.set_xlabel('Time [min]')
        ax3.set_ylabel('Rate')
        fig4, ax4 = plt.subplots()
        ax4.set_xlabel('Time [min]')
        ax4.set_ylabel('Normalized rate')
        ax1.hist(wait_times, bins = 24, edgecolor='black', color = 'gray', alpha = 0.25, density = True, label = 'Single cell data')
        
        
        
        #2# fit with different distribution
        t = np.linspace(0, np.max(wait_times)*1, 10000)
        EMD_dist = dict()
        for dist_name,distr in distribution_list.items():
            if dist_name == 'Weibull':
                P = distr.fit(wait_times,3, floc=0, scale = 30)
            elif dist_name == 'Gamma':
                P = distr.fit(wait_times,floc=0)
            elif dist_name == 'LogNormal':
                P = distr.fit(wait_times,floc=0)
            else:
                P = distr.fit(wait_times)
            wait_times_fit = distr.rvs(*P, size=10000)
            emd_dist = emd_samples(wait_times_fit, wait_times)
            EMD_dist[dist_name] = emd_dist/np.mean(wait_times)
            #TODO = can also implement Hellinger distance
            
            print()
            print("Normalized EMD distance for distribution %s is %.2f%%" % (dist_name,EMD_dist[dist_name]*100))
            print('Optimal parameters for distribution are', P)
            
            
            
            pdf = distr.pdf(t, *P)
            cdf = distr.cdf(t, *P)
            rate = pdf/(1-cdf)
            ax1.plot(t,pdf, lw = 2, label = dist_name)
            ax2.plot(t,cdf, lw = 2, label = dist_name)
            ax3.plot(t,rate, lw = 2, label = dist_name)
            ax4.plot(t,rate/np.max(rate), lw = 2, label = dist_name)
        ax3.legend(prop={'size': 10})
        ax1.legend(prop={'size': 10})
        ax3.set_ylim(0,0.15)
        if wi == 2:
            ax1.set_xlim(25,110)
        plt.show()
    
    """
    print('  Testing our rate functions..')
    t = np.linspace(0.1, np.max(wait_times)*1, 1500)
    plt.figure()
    Gamma_rate_func = np.zeros(t.shape)
    Weibull_rate_func = np.zeros(t.shape)
    Gaussian_rate_func = np.zeros(t.shape)
    Cauchy_rate_func = np.zeros(t.shape)
    Lognormal_rate_func = np.zeros(t.shape)
    for ti,time in enumerate(t):
        Gamma_rate_func[ti] = Gamma_rate(time,10,0.2)
        Weibull_rate_func[ti] = Weibull_rate(time,20,2)
        Gaussian_rate_func[ti] = Gaussian_rate(time, 20, 5)
        Cauchy_rate_func[ti] = Cauchy_rate(time, 20, 1)
        Lognormal_rate_func[ti] = Lognormal_rate(time,2.8,0.2)

    plt.plot(t,Gaussian_rate_func/np.max(Gaussian_rate_func), label = 'Gaussian')  
    plt.plot(t,Weibull_rate_func/np.max(Weibull_rate_func), label = 'Weibull')      
    plt.plot(t,Gamma_rate_func/np.max(Gamma_rate_func), label = 'Gamma')
    plt.plot(t,Cauchy_rate_func/np.max(Cauchy_rate_func), label = 'Cauchy')
    plt.plot(t,Lognormal_rate_func/np.max(Lognormal_rate_func), label = 'Lognormal')
    plt.legend()
    plt.show()
    
    plt.plot(t,Gamma_rate_func)
    plt.show()
    
    plt.plot(t,np.array([Erlang_pdf(ti,20,1) for ti in t]))
    plt.plot(t,np.array([Erlang_cdf(ti,20,1) for ti in t]))
    plt.show()
    """
    
    #https://mathoverflow.net/questions/74545/approximation-to-the-ratio-of-a-gaussian-cdf-to-pdf
    # for the gaussian distribution, rate(t) ~ t for large t
    
    
def calculate_emd(y_sim, alpha, r0):
    """Calculates Earth Movers distance between two distributions
    
    Takes wait-time distribution and calculates the EMD against the theoretical weibull distribution 
    given parameters at the start of the simulation
    
    Args:
        y_sim (list): (n_sim, n_wait_times) list of wait-times for each simulation, 
            for given reaction channel (esc, epi). Obtained from Channels.wait_times attribute
    
    Returns:
        EMD between simulation and true weibull distribution
    """
    y_sim = [item for sublist in y_sim for item in sublist] # flatten list
    y_sim = np.array(y_sim)
    n = 10000
    
    # params for generating theoretical samples (ground truth)
    k = alpha + 1 # shape
    beta = (alpha + 1) * (r0 * gamma((alpha + 2)/(alpha + 1)))**(alpha + 1)
    lam = np.power((alpha+1)/beta, 1/(alpha+1)) # scale
    # generate true samples
    y_true = weibull_min.rvs(k, loc=0, scale=lam, size=n)
    
    d = emd_samples(y_sim, y_true)
    
    return d
    
    
def Weibull_rate(t,mu,alpha):
    #only works for alpha >= 1
    rate = 1/mu
    beta = (alpha) * (rate * gamma((alpha + 1)/(alpha)))**(alpha)
    return np.power(t,alpha-1)*beta

def Erlang_rate(t,k,lbda):
    
    #alternative parametrization with rate: lambda = k/rate
    #only works for k >= 1

    pdf = (lbda**k)*(t**(k-1))*exp(-lbda*t)/factorial(k-1)
    cdf = gammainc(k,lbda*t) #the gammainc function is weird
    rate = pdf/(1-cdf)
    
    return rate

def Gamma_rate(t,alpha,beta):
    #alternative parametrization with rate: beta = alpha/rate
    #only works for alpha >= 1
    pdf = (beta**alpha)*(t**(alpha-1))*exp(-beta*t)/gamma(alpha)
    cdf =  gammainc(alpha,(beta*t)) #the gammainc function is weird
    rate = pdf/(1-cdf)
    return rate

def Gaussian_rate(t, mu, sigma):
    
    #alternative parametrization with rate: mu = 1/rate
    #mp.dps = 50 #set precision of mp function
    if t > mu + 7*sigma: #necessary as the precision of cdf calculation is not precise enough (50 digits precision isnt enough)
        rate = (t-mu)/sigma**2
    else:
        pdf = 1/(sigma*sqrt(2*pi)) * exp(-(1/2) * (t-mu)**2/(sigma**2))
        cdf = 1/2*(1+erf((t-mu)/(sigma*sqrt(2))))
        rate = pdf/(1-cdf)
    return rate

def Cauchy_rate(t,t0,gam):
    #alternative parametrization with rate: t0 = 1/rate
    pdf = 1 / (pi*gam*(1 + ((t-t0)/gam)**2))
    cdf = (1/pi) * atan( (t-t0)/gam ) + 1/2
    rate = pdf/(1-cdf)
    return rate

def Lognormal_rate(t, mu, sigma):
    if t == 0:
        rate = 0
    else:
        pdf = 1/(t*sigma*sqrt(2*pi)) * exp(-(1/2) * (log(t)-mu)**2/(sigma**2))
        cdf = 1/2*(1+erf((log(t)-mu)/(sigma*sqrt(2))))
        rate = pdf/(1-cdf)
    return rate
    
    
def Erlang_pdf(t,k,lbda):
    pdf = (lbda**k)*(t**(k-1))*exp(-lbda*t)/factorial(k-1)
    return pdf

def Erlang_cdf(t,k,lbda):
    cdf = gammainc(k,lbda*t) #the gammainc function is weird
    return cdf

def Gaussian_pdf(t, mu, sigma):
    pdf = 1/(sigma*sqrt(2*pi)) * exp(-(1/2) * (t-mu)**2/(sigma**2))
    return pdf

def Gaussian_cdf(t, mu, sigma):
    cdf = 1/2*(1+erf((t-mu)/(sigma*sqrt(2))))
    return cdf






    

if __name__ == '__main__':
    plt.rcParams.update({'font.size': 16})
    main()