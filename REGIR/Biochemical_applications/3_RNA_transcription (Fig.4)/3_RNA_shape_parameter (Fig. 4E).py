import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from RNA_model import Gillespie_kinetics, get_mRNA_distr_properties

Nsimulations = 200
recompute_mRNA = False
plot_separate = False


def main():
    
    alpha_list = [{'alpha_on':1,'alpha_off':1},{'alpha_on':10,'alpha_off':1},{'alpha_on':1,'alpha_off':10},{'alpha_on':10,'alpha_off':10}]
    #alpha_list = [4]
    
    r_deg = 1
    r_on = 0.8
    r_off = 3
    r_prod = 100
    
    """
    #Nanog:0.0019 0.007 0.80 0.0022
    k_deg = 0.0022
    k_on = 0.0019
    k_off = 0.007
    k_prod = 0.80
    r_on = k_on/k_deg
    r_off = k_off/k_deg
    r_prod = k_prod/k_deg
    """
    
    final_rates = np.array([r_on, r_off, r_prod, r_deg])
    rmin = np.min(final_rates)
    Tend = 7/rmin
    
    if recompute_mRNA:
        mRNA_output = []
        for alpha_params in alpha_list:
            param_list = dict()
            param_list['r_deg'] = 1
            param_list['r_on'] = r_on
            param_list['r_off'] = r_off
            param_list['r_prod'] = r_prod
            param_list['Tend'] = Tend
                
            param_list['alpha_on'] = alpha_params['alpha_on']
            param_list['alpha_off'] = alpha_params['alpha_off']
            param_list['alpha_prod'] = 1
                
            mean_sim, CV, Entropy, Fano, mRNA_population = Gillespie_kinetics(param_list, plot = True, N_simulations = Nsimulations)
            mRNA_output.append(mRNA_population)
            
            # https://thirdorderscientist.org/homoclinic-orbit/2013/10/24/kernel-density-estimation-for-random-variables-with-bounded-support-mdash-the-transformation-trick
            
        np.save('mRNA_output.npy', mRNA_output)
        
    else:
        mRNA_output = np.load('mRNA_output.npy')
        
    plt.rcParams.update({'font.size': 18})
    plt.figure(figsize=(8,5))
    for ai, alpha_params in enumerate(alpha_list):
        data_ = mRNA_output[ai]

        if plot_separate:
            print(alpha_params)
            mean, CV, Entropy, Fano = get_mRNA_distr_properties(data_, pop_is_steady = True, plot = False)
            print('   mean %.2f' % mean)
            print('   CV %.2f' % CV)
            print('   Fano %.2f' % Fano)
            print('   Entropy %.2f' % Entropy)
            print()
                
        #transform data with log function
        data = np.log((data_ + 1))
            
        density = gaussian_kde(data)
        density.covariance_factor = lambda : 0.35 #Smoothing parameter (higher = smoother)
        density._compute_covariance()
        x_vals = np.linspace(np.min(data),np.max(data),200)
        y_vals = density(x_vals)

        x_vals2 = (np.exp(x_vals) - 1)
        y_vals2 = y_vals/(x_vals2 + 1)
            
        #plt.plot(x_vals,y_vals,color = sns.color_palette()[ai+2], lw=2, label = 'alpha = %s' % alpha)
        plt.plot(x_vals2,y_vals2,color = sns.color_palette()[ai+2], lw=3, label = str(alpha_params))
        
        if plot_separate:
            bins = np.arange(0,np.nanmax(data_)+1)
            plt.hist(data_, bins = bins, edgecolor='black', color = 'green', alpha = 0.3, label = 'Simulated', density = True)
            plt.show()
                  
            
    plt.ylabel('PDF')
    plt.xlabel('mRNA counts')
    plt.legend(prop = {'size':14})
    plt.ylim(0,0.05)
    plt.axvline(x=0, color = 'black', lw=2)
    plt.show()
        
        
"""
    t.k.e <- function(x)
    {
        1) log transform the obs
        y = log(x)
        
        2) estimate the density of the pseudo obs with KDE
        k = approxfun(density(y))
        
        3) estimate the density of the original obs
        seq=seq(min(x),max(x),length.out=5e2)
        t.density = as.numeric(vector(length=length(seq)))
        
        for (i in 1:length(seq))
        {
            xx=seq[i]
            t.density[i] = k(log(xx))/xx
        }
        output=list(x=seq,y=t.density)
    }
    
"""
        
        
        
def simBetaPoisson(kon,koff,kprod, num = 1000):
    from scipy.stats import beta, poisson
    rvs = []
    p_array = beta.rvs(kon, koff, size=num)
    for p in p_array:
        n = poisson.rvs(kprod*p, size=1)[0]
        rvs.append(n)
        
    #plt.hist(rvs)
    #plt.show()
    return np.array(rvs)
    
if __name__ == "__main__":
    plt.rcParams.update({'font.size': 16})
    main()