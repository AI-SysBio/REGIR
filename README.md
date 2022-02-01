# A Scalable Gillespie Algorithm for non-Markovian Stochastic Simulations

<img align="right" src="https://raw.githubusercontent.com/Aurelien-Pelissier/REGIR/master/Figures/REGIR.png" width=400>
Discrete stochastic processes are widespread in both nature and human-made systems, with applications across physics, biochemistry, epidemiology, social patterns and finance, just to name a few. In the majority of these systems, the dynamics cannot be properly described with memoryless (or Markovian) interactions, and thus require the use of non-Markovian simulations. This repository contains an implementattion of a general and scalable framework to simulate non-Markovian stochastic systems with arbitrary inter-event time distribution and accuracy. The algorithm is referred to as the Rejection Gillespie algorithm for non-Markovian Reactions (REGIR) [1].

&nbsp;


The algorithm efficiency and accuracy lies in its rejection approach, whereby processes that include undesired reactions are computed and subsequently rejected. See the corresponding paper for details. With the current implementation, the following distribution are available.  

      Exponential distribution:
      
      Normal distrubution:
      
      LogNormal distribution:
      
      Gamma distribution:
      
      Weibull distribution:
      
      Cauchy distribution:
      
      

      
Note that monotolically decreasing distribution, such as Weibull ($\alpha < 1$), gamma ($\alpha < 1$) or power laws, are not available in the current implementation of this repository, as theses can be more elegantly and efficiently simulated with the Laplace Gillespie algorithm [2]. Feel free to drop me an email if you would be interrested in adding the Laplace Gillespie to this repository, or adding a new distribution of your interrest.
        
        
### Running the code
To launch the simulation, run `main.py`. Running the program requires python3 with its standard libraries such as numpy or matplotlib. The program plot the simulation results from the model, along with the experimental data from litterature, available in the folder `Exp_data/`. One Gillespie GC run should take less than 1min.


## References

[1] Pélissier, A, Phan, M, et al. Practical and scalable simulations of non-Markovian stochastic processes. Proceedings of the National Academy of Sciences (2022)

[2] Masuda, Naoki, and Luis EC Rocha. "A Gillespie algorithm for non-Markovian stochastic processes." Siam Review 60.1 (2018): 95-115.
