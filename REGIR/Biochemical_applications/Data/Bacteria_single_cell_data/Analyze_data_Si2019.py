import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os,sys

from scipy.stats import lognorm, expon, gengamma, weibull_min, norm, cauchy

def main():
    B_subtilis_data_file = 'mbio.02205-19-data.xlsx'
    """
    B. subtilis single-cell data for both growth and the cell cycle.
    """
    
    B_C_data_file = 'Si2019.xlsx'
    """ 
    provide a complete single-cell growth and cell cycle data 
    in both E. coli and B. subtilis, tracking not only the division cycles 
    but also the replication cycles. 
    
    We implemented the Helmstetter-Cooper mode
    that is often interpreted to mean that initiation triggers 
    division after a fixed elapsed time τcyc = C + D
    """
    
    alpha_list = []
    mean_list = []
    for n_sheet in range(1,8):
        BC_data = pd.read_excel(B_C_data_file, sheet_name=n_sheet)
        #print(B_data)
        print(BC_data.columns)
        
        #birth_time = B_data['birth time (minute)'].to_numpy()
        #division_time = B_data['division time (minute)'].to_numpy()
        #initiation_time = B_data['initiation time (minute)'].to_numpy()
        #growth_rate = B_data['growth rate (1/hour)'].to_numpy()
        #generation_time = division_time - birth_time
        generation_time = BC_data['generation time (minute)'].to_numpy()
        B_period = BC_data['B period (minute)'].to_numpy()
        C_period = BC_data['C period (minute)'].to_numpy()
        D_period = BC_data['D period (minute)'].to_numpy()
        
        
        data = generation_time
        #data = D_period
        
        P = weibull_min.fit(data,3, floc=0, scale = 30)
        alpha_list.append(P[0]-1)
        mean_list.append(np.mean(data))
        P2 = cauchy.fit(data)
        print(P)
        print('(alpha+1, x0, scale)')
        x = np.linspace(0, 1.1*np.max(data), 300)
        pdf = weibull_min.pdf(x, *P)
        pdf2 = cauchy.pdf(x, *P2)
        
        plt.hist(data, density = True, color = 'black', bins = 24, alpha=0.5)
        plt.plot(x, pdf2, 'g',lw=2, label='normal') 
        plt.plot(x, pdf, 'r',lw=2, label='Weibull') 
        plt.legend()
        plt.ylabel('PDF')
        plt.xlabel('Time to division [min]')
        plt.show()
    
        """
        plt.hist(growth_rate, density = True, color = 'black', bins = 30)
        plt.ylabel('PDF')
        plt.xlabel('Growthrate [1/h]')
        plt.show()
        """
        
    print(alpha_list)
    print(mean_list)

    plt.figure(figsize = (1.5,5))
    plt.scatter(np.zeros(7), alpha_list,color='black')
    plt.ylim(0,7)
    plt.xticks([])
    plt.ylabel('alpha')
    plt.axhline(y=np.median(alpha_list), color = 'black')
    plt.show()
    
    plt.scatter(mean_list,alpha_list, color = 'black')
    plt.xlabel('mean time to division')
    plt.ylabel('alpha parameter')
    plt.show()



if __name__ == '__main__':
    plt.rcParams.update({'font.size': 15})
    main()