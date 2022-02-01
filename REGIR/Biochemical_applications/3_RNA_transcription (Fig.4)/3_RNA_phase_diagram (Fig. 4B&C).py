import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import random as rd
import os

from RNA_model import Gillespie_kinetics, get_mRNA_distr_properties


nres = 5
recompute_phase_diagram = True
N_simulations = 20
rd.seed(401)
    
    
def main():
    
    if not os.path.isdir("Phase_diagrams"):
        os.makedirs("Phase_diagrams")
    

    #1-phase diagram
    phase_diagram_betapoisson(nres = nres) #makovian only
    
    alpha_on = 10
    alpha_off = 10
    alpha_prod = 1
    phase_diagram(alpha_off = alpha_off, alpha_on = alpha_on, alpha_prod = alpha_prod, nres = nres, plot_distrib = True, log_space = False)
    
    

def phase_diagram_betapoisson(nres = 4, plot_distrib = True, log_space = False):
    
    """
    [ref] Vu, Trung Nghia, et al. "Beta-Poisson model for single-cell 
    RNA-seq data analyses." Bioinformatics 32.14 (2016): 2128-2135.
    """
    
    
    r_prod_bounds = [100,100]
    r_on_bounds = [1e-3,1e1]
    r_off_bounds = [1e-3,1e1]
    
    #for entropy plot
    r_on_bounds = [0.2,5]
    r_off_bounds = [0.2,5]
    
    if r_prod_bounds[1] > r_prod_bounds[0]:
        if log_space:
            r_prod_array = np.logspace(np.log10(r_prod_bounds[0]),np.log10(r_prod_bounds[1]), num = nres)
        else:
            r_prod_array = np.linspace(r_prod_bounds[0],r_prod_bounds[1], num = nres)
    else:
        r_prod_array = [r_prod_bounds[0]]
        
    if log_space:
        r_on_array = np.logspace(np.log10(r_on_bounds[0]),np.log10(r_on_bounds[1]), num = nres)
        r_off_array = np.logspace(np.log10(r_off_bounds[0]),np.log10(r_off_bounds[1]), num = nres)
    else:
        r_on_array = np.linspace(r_on_bounds[0],r_on_bounds[1], num = nres)
        r_off_array = np.linspace(r_off_bounds[0],r_off_bounds[1], num = nres)
        
    if recompute_phase_diagram:
    
        CV_arr = np.zeros((len(r_prod_array),len(r_off_array),len(r_on_array)))
        Entropy_arr = np.zeros((len(r_prod_array),len(r_off_array),len(r_on_array)))
        Fano_arr = np.zeros((len(r_prod_array),len(r_off_array),len(r_on_array)))
        n_sim_tot = len(r_prod_array) * len(r_on_array) * len(r_off_array)
        n_sim = 0
        for pi, r_prod in enumerate(r_prod_array):
            for oni, r_on in enumerate(r_on_array):
                for offi, r_off in enumerate(r_off_array):
                    n_sim += 1
                    print('  Simulation %s/%s ...' % (n_sim,n_sim_tot))
                    param_list = dict()
                    r_deg = 1
                    param_list['r_deg'] = r_deg
                    param_list['r_on'] = r_on
                    param_list['r_off'] = r_off
                    param_list['r_prod'] = r_prod
                    print('   ',param_list)
                    
                    rmin = np.min([r_on,r_off,r_prod,r_deg])
                    Tend = 5/rmin # important to reach steady state and correctly assess distribution
                    param_list['Tend'] = Tend
                    
                    mRNA_count = simBetaPoisson(r_on,r_off,r_prod, num = 10000)
                    mean, CV, Entropy, Fano = get_mRNA_distr_properties(mRNA_count, plot = plot_distrib, pop_is_steady = True)
                    
                    CV_arr[pi,offi,oni] = CV
                    Fano_arr[pi,offi,oni] = Fano
                    Entropy_arr[pi,offi,oni] = Entropy
        
                    
        np.save('Phase_diagrams/CV_arr_bp_%s.npy' % nres,CV_arr)
        np.save('Phase_diagrams/Fano_arr_bp_%s.npy' % nres,Fano_arr)
        np.save('Phase_diagrams/Entropy_arr_bp_%s.npy' % nres,Entropy_arr)
    
    else:
        
        CV_arr = np.load('Phase_diagrams/CV_arr_bp_%s.npy' % nres)
        Fano_arr = np.load('Phase_diagrams/Fano_arr_bp_%s.npy' % nres)
        Entropy_arr = np.load('Phase_diagrams/Entropy_arr_bp_%s.npy' % nres)
    
    
    plot_results(CV_arr,Fano_arr,Entropy_arr,r_on_array,r_off_array,log_space_xy = log_space)
                
    

        
        
def phase_diagram(alpha_off = 1, alpha_on = 1, alpha_prod = 1, nres = 4, plot_distrib = True, log_space = False):
    
    """
    We find that the results of our exact solution can be interpreted in 
    the most natural and transparent way when we measure time in units of 1/kd,
    i.e., in terms of the mRNA lifetime. Thus we will use the three dimensionless 
    ratios kb/kd, cf/kd, and cb/kd to organize our results.
    
    We note that the mean number of mRNA in the steady state is given by the product 
    of cf / (cf + cb), the fraction of the time the gene is in the activated state,
    and kb /kd, the mean value of mRNA if the gene is “on.”
    
    - The ratio kb /kd sets the scale for the number of mRNA and increasing it 
      extends the range over which P(m) is appreciable without a significant change of shape. 
      
    - The remaining ratios cf /kd and cb /kd determine the shape of the distribution.
    """
    r_prod_bounds = [100,100]
    r_on_bounds = [1e-3,1e1]
    r_off_bounds = [1e-3,1e1]
    
    #for entropy plot
    r_on_bounds = [0.2,5]
    r_off_bounds = [0.2,5]
    
    #for CV plot
    #log_space = True
    #r_on_bounds = [1e-2,1e1]
    #r_off_bounds = [1e-2,1e1]
    
    if r_prod_bounds[1] > r_prod_bounds[0]:
        if log_space:
            r_prod_array = np.logspace(np.log10(r_prod_bounds[0]),np.log10(r_prod_bounds[1]), num = nres)
        else:
            r_prod_array = np.linspace(r_prod_bounds[0],r_prod_bounds[1], num = nres)
    else:
        r_prod_array = [r_prod_bounds[0]]
        
    if log_space:
        r_on_array = np.logspace(np.log10(r_on_bounds[0]),np.log10(r_on_bounds[1]), num = nres)
        r_off_array = np.logspace(np.log10(r_off_bounds[0]),np.log10(r_off_bounds[1]), num = nres)
    else:
        r_on_array = np.linspace(r_on_bounds[0],r_on_bounds[1], num = nres)
        r_off_array = np.linspace(r_off_bounds[0],r_off_bounds[1], num = nres)
    
    if recompute_phase_diagram:
        CV_arr = np.zeros((len(r_prod_array),len(r_off_array),len(r_on_array)))
        Entropy_arr = np.zeros((len(r_prod_array),len(r_off_array),len(r_on_array)))
        Fano_arr = np.zeros((len(r_prod_array),len(r_off_array),len(r_on_array)))
        n_sim_tot = len(r_prod_array) * len(r_on_array) * len(r_off_array)
        n_sim = 0
        for pi, r_prod in enumerate(r_prod_array):
            for oni, r_on in enumerate(r_on_array):
                for offi, r_off in enumerate(r_off_array):
                    n_sim += 1
                    print('  Simulation %s/%s ...' % (n_sim,n_sim_tot))
                    param_list = dict()
                    r_deg = 1
                    param_list['r_deg'] = r_deg
                    param_list['r_on'] = r_on
                    param_list['r_off'] = r_off
                    param_list['r_prod'] = r_prod
                    param_list['alpha_on'] = alpha_on
                    param_list['alpha_off'] = alpha_off
                    param_list['alpha_prod'] = alpha_prod
                    print('   ',param_list)
                    
                    rmin = np.min([r_on,r_off,r_prod,r_deg])
                    Tend = 5/rmin # important to reach steady state and correctly assess distribution
                    param_list['Tend'] = Tend
                    mean, CV, Entropy, Fano, mRNApop_steady = Gillespie_kinetics(param_list, plot = False, plot_distrib = plot_distrib, N_simulations = N_simulations)
                    
                    CV_arr[pi,offi,oni] = CV
                    Fano_arr[pi,offi,oni] = Fano
                    Entropy_arr[pi,offi,oni] = Entropy
                
        np.save('Phase_diagrams/CV_arr_%s_%s_%s.npy' % (alpha_on, alpha_off, nres),CV_arr)
        np.save('Phase_diagrams/Fano_arr_%s_%s_%s.npy' % (alpha_on, alpha_off, nres),Fano_arr)
        np.save('Phase_diagrams/Entropy_arr_%s_%s_%s.npy' % (alpha_on, alpha_off, nres),Entropy_arr)
    
    else:
        
        CV_arr = np.load('Phase_diagrams/CV_arr_%s_%s_%s.npy' % (alpha_on, alpha_off, nres))
        Fano_arr = np.load('Phase_diagrams/Fano_arr_%s_%s_%s.npy' % (alpha_on, alpha_off, nres))
        Entropy_arr = np.load('Phase_diagrams/Entropy_arr_%s_%s_%s.npy' % (alpha_on, alpha_off, nres))    
    
    plot_results(CV_arr,Fano_arr,Entropy_arr,r_on_array,r_off_array,log_space_xy = log_space)
                
                

def plot_results(CV_arr,Fano_arr,Entropy_arr,r_on_array,r_off_array,log_space_xy = False, log_space_c = False):      
    
    mean_arr = Fano_arr/np.power(CV_arr,2)
    y = r_on_array
    x = r_off_array

    
    if log_space_xy:
        xnew = np.logspace(np.log10(np.min(r_off_array)), np.log10(np.max(r_off_array)), 100)
        ynew = np.logspace(np.log10(np.min(r_on_array)), np.log10(np.max(r_on_array)), 100)
    else:
        xnew = np.linspace(np.min(r_off_array), np.max(r_off_array), 100)
        ynew = np.linspace(np.min(r_on_array), np.max(r_on_array), 100)
    ticks_pos_on = np.linspace(0,len(xnew)-1, num = 5).astype(int)
    ticks_pos_off = np.linspace(0,len(ynew)-1, num = 5).astype(int)
    
    
    """ mean"""
    log_space_c = False
    if log_space_c:
        f = interpolate.interp2d(x, y, np.log(mean_arr[0,:,:]), kind='cubic')
        levels = np.linspace(np.min(np.log(mean_arr[0,:,:])),np.max(np.log(mean_arr[0,:,:])),num = 15)
    else:
        f = interpolate.interp2d(x, y, mean_arr[0,:,:], kind='cubic')
        levels = np.linspace(np.min(mean_arr[0,:,:]),np.max(mean_arr[0,:,:]),num = 15)
        levels = np.linspace(3,100,num = 15)
    znew = f(xnew, ynew)
    
    plt.figure(figsize=(8.5,7))
    if log_space_c:
        plt.imshow(np.transpose(znew), cmap = 'inferno', origin = 'lower')
    else:
        plt.imshow(np.transpose(znew), cmap = 'inferno', origin = 'lower', vmin = 3, vmax=100)
    cbar = plt.colorbar()
    cbar.set_label('Mean',size=18)
    tick_values = np.linspace(np.min(znew),np.max(znew), num = 5)
    cbar.set_ticks(tick_values)
    if log_space_c:
        cbar.set_ticklabels(np.round(np.exp(tick_values),decimals=2))
    else:
        tick_values = np.linspace(3,100, num = 5)
        cbar.set_ticks(tick_values)
        cbar.set_ticklabels(np.round(tick_values,decimals=2))
    plt.contour(np.transpose(znew), levels, colors = 'black')
    plt.yticks(ticks = ticks_pos_on, labels = np.round(ynew[ticks_pos_on],decimals = 2))
    plt.xticks(ticks = ticks_pos_off, labels = np.round(xnew[ticks_pos_off],decimals = 2))
    plt.ylabel(r'$\lambda_{on} \ / \ \lambda_{deg}$', size = 19)
    plt.xlabel(r'$\lambda_{off} \ / \ \lambda_{deg}$', size = 19)
    plt.show()
    
    
    
    """ CV """
    log_space_c = True
    if log_space_c:
        f = interpolate.interp2d(x, y, np.log(CV_arr[0,:,:]), kind='cubic')
        levels = np.linspace(np.min(np.log(CV_arr[0,:,:])),np.max(np.log(CV_arr[0,:,:])),num = 15)
        levels = np.linspace(np.log(0.1),np.log(3),num = 15)
    else:
        f = interpolate.interp2d(x, y, CV_arr[0,:,:], kind='cubic')
        levels = np.linspace(np.min(CV_arr[0,:,:]),np.max(CV_arr[0,:,:]),num = 15)
    znew = f(xnew, ynew)
    
    plt.figure(figsize=(8.5,7))
    if log_space_c:
        plt.imshow(np.transpose(znew), cmap = 'inferno', origin = 'lower', vmin = np.log(0.1), vmax = np.log(3))
    else:
        plt.imshow(np.transpose(znew), cmap = 'inferno', origin = 'lower')
    cbar = plt.colorbar()
    cbar.set_label('CV',size=18)
    tick_values = np.linspace(np.min(znew),np.max(znew), num = 5)
    cbar.set_ticks(tick_values)
    if log_space_c:
        tick_values = np.linspace(np.log(0.1),np.log(3), num = 5)
        cbar.set_ticks(tick_values)
        cbar.set_ticklabels(np.round(np.exp(tick_values),decimals=2))
    else:
        cbar.set_ticklabels(np.round(tick_values,decimals=2))
    plt.contour(np.transpose(znew), levels, colors = 'black', linestyles = 'solid')
    plt.yticks(ticks = ticks_pos_on, labels = np.round(ynew[ticks_pos_on],decimals = 2))
    plt.xticks(ticks = ticks_pos_off, labels = np.round(xnew[ticks_pos_off],decimals = 2))
    plt.ylabel(r'$\lambda_{on} \ / \ \lambda_{deg}$', size = 19)
    plt.xlabel(r'$\lambda_{off} \ / \ \lambda_{deg}$', size = 19)
    plt.show()
    
    
    """ Fano factor """
    log_space_c = True
    if log_space_c:
        f = interpolate.interp2d(x, y, np.log(Fano_arr[0,:,:]), kind='cubic')
    else:
        f = interpolate.interp2d(x, y, Fano_arr[0,:,:], kind='cubic')
    znew = f(xnew, ynew)
    
    plt.figure(figsize=(8.5,7))
    plt.imshow(np.transpose(znew), cmap = 'inferno', origin = 'lower')
    cbar = plt.colorbar()
    cbar.set_label('Fano factor',size=18)
    tick_values = np.linspace(np.min(znew),np.max(znew), num = 5)
    cbar.set_ticks(tick_values)
    if log_space_c:
        cbar.set_ticklabels(np.round(np.exp(tick_values),decimals=2))
    else:
        cbar.set_ticklabels(np.round(tick_values,decimals=2))
    plt.contour(np.transpose(znew), colors = 'black')
    plt.yticks(ticks = ticks_pos_on, labels = np.round(ynew[ticks_pos_on],decimals = 2))
    plt.xticks(ticks = ticks_pos_off, labels = np.round(xnew[ticks_pos_off],decimals = 2))
    plt.ylabel(r'$\lambda_{on} \ / \ \lambda_{deg}$', size = 19)
    plt.xlabel(r'$\lambda_{off} \ / \ \lambda_{deg}$', size = 19)
    plt.show()
    
    
    
    """ Entropy """
    #obviously logspace_c = False
    z = Entropy_arr[0,:,:]
    zs = z
    #zs = sgolay2.SGolayFilter2(window_size=17, poly_order=2)(z) #for nres = 30
    #zs = sgolay2.SGolayFilter2(window_size=3, poly_order=1)(z)#for nres = 5
    #zs = sgolay2.SGolayFilter2(window_size=3, poly_order=1)(z)#for nres = 10

    f = interpolate.interp2d(x, y, zs, kind='cubic')
    znew = f(xnew, ynew)
    
    levels = np.arange(0,4,step=0.1)
    
    plt.figure(figsize=(8.5,7))
    if log_space_xy:
        plt.imshow(np.transpose(znew), cmap = 'inferno', origin = 'lower', vmax = 3)
    else:
        plt.imshow(np.transpose(znew), cmap = 'inferno', origin = 'lower', vmax = 3, vmin = 2)
    cbar = plt.colorbar()
    cbar.set_label('Entropy',size=18)
    if log_space_xy:
        tick_values = np.linspace(np.min(Entropy_arr),3, num = 5)
    else:
        tick_values = np.linspace(2,3, num = 5)
    cbar.set_ticks(tick_values)
    cbar.set_ticklabels(np.round(tick_values,decimals = 2))
    plt.contour(np.transpose(znew), levels, colors = 'black')
    plt.yticks(ticks = ticks_pos_on, labels = np.round(ynew[ticks_pos_on],decimals = 2))
    plt.xticks(ticks = ticks_pos_off, labels = np.round(xnew[ticks_pos_off],decimals = 2))
    plt.ylabel(r'$\lambda_{on} \ / \ \lambda_{deg}$', size = 19)
    plt.xlabel(r'$\lambda_{off} \ / \ \lambda_{deg}$', size = 19)
    plt.show()
    
        
        
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
    plt.rcParams.update({'font.size': 20})
    main()