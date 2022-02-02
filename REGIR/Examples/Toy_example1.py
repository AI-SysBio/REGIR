import REGIR as gil

#Set the simulation parameters:
class param:
	Tend = 10		#Length of the simulation
	unit = 'h'		#Unit of time (is used for plotting purpose only)
	N_simulations = 20	#The simulation results should be averaged over many trials
	timepoints = 100	#Number of timepoints to record (make surethat this number isnt too big)

r1 = 1
r2 = 4
r3 = 0.03
alpha1 = 20
alpha2 = 5
  
#Define the reaction chanels:
reaction1 = gil.Reaction_channel(param,rate=r1, shape_param=alpha1, distribution = 'Gamma', name = 'A -> B')
reaction1.reactants = ['A']
reaction1.products = ['B']	
reaction2 = gil.Reaction_channel(param,rate=r2, shape_param=alpha2, distribution = 'Weibull', name = 'A -> A+C')
reaction2.reactants = ['B']
reaction2.products = ['A','C']	
reaction3 = gil.Reaction_channel(param,rate=r3, name = 'A + B -> 0')
reaction3.reactants = ['A','B']
reaction3.products = []
	

#Define the initial population of reactants:
N_init = dict()
N_init['A'] = 300
N_init['B'] = 0
N_init['C'] = 0

#Initialize and run the Gillepsie simulation:
reaction_channel_list = [reaction1, reaction2, reaction3]
G_simul = gil.Gillespie_simulation(N_init,param)
G_simul.reaction_channel_list = reaction_channel_list
populations = G_simul.run_simulations(param.Tend)
G_simul.plot_inter_event_time_distribution()
G_simul.plot_populations()