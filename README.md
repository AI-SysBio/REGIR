# A Scalable Gillespie Algorithm for non-Markovian Stochastic Simulations

<img align="right" src="https://raw.githubusercontent.com/Aurelien-Pelissier/REGIR/master/Figures/REGIR.png" width=400>
Discrete stochastic processes are widespread in both nature and human-made systems, with applications across physics, biochemistry, epidemiology, social patterns and finance, just to name a few. In the simplest case, these processes are memoryless (or Markovian), with future occurrences predictable based solely on the present state of the system and with exponentially distributed interevent times. However, stochastic systems describing majority of applications are empirically known to exhibit properties of memory, an inherently non-Markovian feature. This repository contains an implementattion of a general and scalable framework to simulate non-Markovian stochastic systems with arbitrary inter-event time distribution and accuracy. The algorithm is referred to as the Rejection Gillespie algorithm for non-Markovian Reactions (REGIR) [1].

&nbsp;


The model is based on a system of 11 stochastic interactions, implemented with a modified Gillepsie algorithm that can account for the individual properties of each agents [2].

      8 types of reactants are considered:
      - Centroblast (Dark Zone) = CB
      - Centrocytes (Light Zone) = CC
      - Selected Centrocytes (Light Zone) = CCsel
      - Bound Centrocytes (Light Zone) = [CCTC]
      - Free T follicular helper (Light Zone) = Tfh
        (Plus 3 additional cell types, leaving the GC)
      - Memory cells (Outside GC) = MBC
      - Plasma cells (Outside GC) = PC
      - Dead cells = 0 
      
      11 reactions are considered:
      - Cell entering the GC:        0 -> CB
      - Centrocyte apoptosis:        CC -> 0
      - Centroblast migration:       CB -> CC
      - Centrocyte unbinding:        [CCTC] -> CC + TC
      - Centrocyte recirculation:    CC -> CB
      - Centrocyte exit:             CC -> MBC or PC
      - Centrocyte Tfh binding:      CC + TC = [CCTC]
      - Tfh switch:                  [CC1TC] + CC2 -> CC1 + [CC2TC]
      - Centroblast division:        CB -> 2CB
      - Centroblast apoptosis:       CB -> 0
      - Centrocyte antigen uptake:   CC -> CC
        
        
### Running the code
To launch the simulation, run `main.py`. Running the program requires python3 with its standard libraries such as numpy or matplotlib. The program plot the simulation results from the model, along with the experimental data from litterature, available in the folder `Exp_data/`. One Gillespie GC run should take less than 1min.


## References

[1] Pélissier, A, Phan, M, et al. Practical and scalable simulations of non-Markovian stochastic processes. Proceedings of the National Academy of Sciences (2022)
