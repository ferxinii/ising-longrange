# Monte Carlo simulation of the Ising model with super long-range interactions ⚡️
This was initally developed as an optional project for my Statistical Physics course. However, I revisited the project later and improved it. This repository includes the Python code that runs the simulations, as well as the MATLAB code used in the numerical analysis of the analytical model.

In **numeric_matlab** you will find everything related to the numerical exploration, including the plots used for the project.

In **results** you will find the results of running **main.py**, as well as the results of previous runs. Additionally, a MATLAB script reads those previous runs to generate the plot that is used for the project.

## Running the code
1. simply change the desired parameters in the main function. These include the number of atoms and number of Monte Carlo steps, as well as the temperature range. Don't go overboard or it might be too slow.
2. Then you can just run it. An indicator lets you know how many steps of temperature you have done so you can keep the progress.
3. Once finished, you will see a plot with the results of the run. This same plot and the data is automatically saved inside **simulation_results**.
