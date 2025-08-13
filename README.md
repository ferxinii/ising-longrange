# Ising model with super long-range interactions ⚡️
We consider an Ising model with hamiltonian:
```math
\mathcal{H}(\{ s_i\}) = \frac{-J}{N} \sum_{i=1}^{N-1} \sum_{j>i}^N s_i s_j - H \sum_{i=1}^Ns_i
```

Here, $$J$$ is the strenght of interaction, and $$H$$ an external magnetic field.

This repository includes:
-  pdf with an **analytic analysis** of the model and an exploration of the **Monte Carlo simulation** results.
-  Python code to run MC simulations of the model (for the case of zero field).

This was initally developed as a project for my Statistical Physics course, with a later revision to improve it.

## Running the code
1. Simply run *main.py*. Input the number of spins and the number of Monte Carlo steps to perform. I recommend starting with 50 spins and 10000 steps.
2. Once finished, a folder is created with the results of the run: A file with the mean magnetization per spin and total energy of each of the MC samples, for each of the temperatures. Be careful with over-writing previous results.

An example of the output is found in the folder N50.


## Processed results from MC simulations
<p align="center">
<img src="./images/summary.png" alt="Summary of MC results" width="400" height="auto" />
<img src="./images/m_evolution.jpg" alt="Moving average of m" width="400" height="auto">
</p>


## Analytical results (see pdf)
Some of the relevant equations found in the project:

 ``` math
\begin{matrix}
E = \frac{J}{2} - HM - \frac{J}{2N} M^2  &  m = \tanh{ \left( \frac{J}{k_BT} (m + \frac{H}{J}) \right)} \\ 
c_V(T,H) = \frac{J}{T} \frac{1-m^2}{m^2 - \left( 1-\frac{k_B T}{J} \right)} \left( m^2 + \frac{H}{J}m\right)  &  \chi _M(T, H) = \frac{1}{J} \frac{1-m^2}{m^2-\left( 1- \frac{k_B T}{J}\right)}
\end{matrix}
```


<p align="center">
<img src="./images/surfaces_TH.png" alt="Surfaces as a function of T and H" width="400" height="auto"/>
</p>
