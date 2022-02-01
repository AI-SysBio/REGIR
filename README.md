# A Scalable Gillespie Algorithm for non-Markovian Stochastic Simulations

<img align="right" src="https://raw.githubusercontent.com/Aurelien-Pelissier/REGIR/master/Figures/REGIR.png" width=400>
Discrete stochastic processes are widespread in both nature and human-made systems, with applications across physics, biochemistry, epidemiology, social patterns and finance, just to name a few. In the majority of these systems, the dynamics cannot be properly described with memoryless (or Markovian) interactions, and thus require the use of numerical tools for simulating these non-Markovian dynamics. This repository contains an implementattion of a general and scalable framework to simulate non-Markovian stochastic systems with arbitrary inter-event time distribution and accuracy. The algorithm is referred to as the Rejection Gillespie algorithm for non-Markovian Reactions (REGIR) [1].

&nbsp;


With the current implementation, each distributions are characterised by their rate and a shape parameter as follow:

      Exponential:
          - rate: 1/mean
          - shape parameter: None
      
      Normal:
          - rate: 1/mean
          - shape: std/mean
      
      LogNormal:
          - rate: 1/mean
          - shape: std/mean
          
      Gamma:
          - rate: 1/mean
          - shape: α >= 1 (https://en.wikipedia.org/wiki/Gamma_distribution)
          
      Weibull:
          - rate: 1/mean
          - shape: k >= 1 (https://en.wikipedia.org/wiki/Weibull_distribution)
          
      Cauchy:
          - rate: 1/median
          - shape: γ (https://en.wikipedia.org/wiki/Cauchy_distribution)
      

      
Note that monotolically decreasing distributions, such as Weibull (k < 1), gamma (α < 1) or power laws, are not available in the current implementation of this repository, as theses can be more elegantly and efficiently simulated with the Laplace Gillespie algorithm [2]. Feel free to drop me an email if you would be interrested in adding the Laplace Gillespie or any other distributions of your interrest to this repository.
        
        
### Simulating a non-Markovian system

To launch a Gillespie simulation, we do as follow:

      GAGAG
          - rate: 1/mean
          - shape parameter: None
      
      Normal:
          - rate: 1/mean
          - shape: std/mean

Keep in mind that non-Markovian simulations are only available for reaction chanels with a single reactant, as the definition of inter-event time distribution is ambigious for channels with multiple reactants. If a chanel is defined with 0 or more than one reactant, it will be considered as a Poisson process. Simple examples, as well as the three biochemical systems described in the paper (Cell division, differentiation and RNA transcription) are provided in the example folder. 


## References

[1] Pélissier, A, Phan, M, et al. Practical and scalable simulations of non-Markovian stochastic processes. Proceedings of the National Academy of Sciences (2022)

[2] Masuda, Naoki, and Luis EC Rocha. "A Gillespie algorithm for non-Markovian stochastic processes." Siam Review 60.1 (2018): 95-115.
