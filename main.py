"""This code performs a Monte Carlo simulation of the Ising model considering long-range interactions for a range of T.
Improved version of a code from June 2022, rewritten in January 2024. """

import copy
import numpy as np
import matplotlib.pyplot as plt
import os
import numba as nb  # This improves performance, a lot


# -----------------------------------------------------------------------------------------------------------
# ----------------- CLASS DEFINITIONS -----------------------------------------------------------------------
class State:
    """Characterized by the number of particles, their spin (np_array), and the hamiltonian.
    It initializes by default to a random state. Compiled with numba for performance."""
    N: int
    spins: np.ndarray
    ham: float

    def __init__(self, N: int):
        self.N = N
        self.random_state()  # Updates spins and ham

    def random_state(self):
        """Update "spins" and "ham" with a random state."""
        self.spins = np.random.choice([1, -1], size=self.N)
        self.ham = self.hamiltonian()

    def change_spin(self):
        """Changes state of a random spin, and updates ham."""
        rand_id = np.random.randint(0, self.N)
        self.spins[rand_id] *= -1
        self.ham = self.hamiltonian()

    def hamiltonian(self) -> float:
        """Computes the hamiltonian."""
        return hamiltonian_optimized(self.N, self.spins)


class Ensemble:
    """An Ensemble is a collection of M states (each of N particles) used for sampling. Contains the number of samples
    already stored, an array of arrays with the spins of each sample, and the hamiltonian of each sample."""
    M: int
    spins: np.ndarray
    ham: np.ndarray
    ii: int

    def __init__(self, M: int, N: int):
        """M is the number of MC steps, N is the number of particles. Initialized to 0s."""
        self.M = M
        self.spins = np.empty((M, N))
        self.ham = np.empty(M)
        self.ii = 0  # This is to know what is the last spins ii have stored, counter of calls to add_ensemble.

    def add_sample(self, state: State):
        """Adds into the ensemble the desired state, which is an object of the class State."""
        if self.ii < self.M:
            self.spins[self.ii][:] = state.spins
            self.ham[self.ii] = state.ham
            self.ii += 1
        else:
            print("Already filled the ensemble!")
            exit(1)


class Statistics:
    """Class to store and compute the statistics of a given Ensemble, for a range of T."""
    def __init__(self, Ti: float, Tf: float, it_temp: int, M: int, N: int):
        self.M = M
        self.N = N
        self.ii = 0  # Counter of the calls made to average_ensemble
        self.it_temp = it_temp
        self.T_vec = np.linspace(Ti, Tf, it_temp)
        self.mag = np.zeros(it_temp)  # Mean magnetic moment per spin, at each T
        self.MAG = np.zeros(it_temp)  # Mean total magnetization, at each T
        self.MAG2 = np.zeros(it_temp)  # Mean of the squared total magnetization, at each T
        self.E = np.zeros(it_temp)  # Mean total energy, at each T
        self.E2 = np.zeros(it_temp)  # Mean of the square of the total energy, at each T
        self.cv = np.zeros(it_temp)
        self.xm = np.zeros(it_temp)

    def average_ensemble(self, ensemble: Ensemble):
        """Calculates statistical averages from the given Ensemble."""
        if self.ii < self.it_temp:
            MAG_aux = np.abs(np.sum(ensemble.spins, axis=1))  # total magnetization of each sample
            self.MAG[self.ii] = 1 / self.M * np.sum(MAG_aux)
            self.MAG2[self.ii] = 1 / self.M * np.sum(np.square(MAG_aux))
            self.mag[self.ii] = 1 / self.N * self.MAG[self.ii]
            self.E[self.ii] = 1 / self.M * np.sum(ensemble.ham)
            self.E2[self.ii] = 1 / self.M * np.sum(np.square(ensemble.ham))
            self.ii += 1
        else:
            print("Already filled the statistics!")
            exit(1)

    def compute_fluctuations(self):
        """Should only be called once all the data has been stored, for all temperatures."""
        if self.ii == self.it_temp:
            beta = 1 / self.T_vec
            self.cv = 1 / self.N * np.square(beta) * (self.E2 - np.square(self.E))
            self.xm = 1 / self.N * beta * (self.MAG2 - np.square(self.MAG))
        else:
            print("Could not compute fluctuations. Only call when all temperatures have been stored.")
            exit(1)

    def save_file(self):
        """Save the final data in a document."""
        save_data = open("./simulation_results/data_N"+str(self.N), "w")
        save_data.write(' '.join([str(element) for element in self.T_vec]))
        save_data.write('\n')
        save_data.write(' '.join([str(element) for element in self.mag]))
        save_data.write('\n')
        save_data.write(' '.join([str(element) for element in self.E]))
        save_data.write('\n')
        save_data.write(' '.join([str(element) for element in self.cv]))
        save_data.write('\n')
        save_data.write(' '.join([str(element) for element in self.xm]))
        save_data.write('\n')
        text = ['N: ', str(self.N), '   it_mc: ', str(self.M), '   T0: ', str(self.T_vec[0]),
                '    Tf: ', str(self.T_vec[-1]), '   it_temp: ', str(self.it_temp)]
        save_data.writelines(text)
        save_data.close()


# -----------------------------------------------------------------------------------------------------------
# ----------------- FUNCTION DEFINITIONS --------------------------------------------------------------------
@nb.njit
def hamiltonian_optimized(N: int, spins: np.ndarray) -> float:
    """Returns hamiltonian of a given set of spins. Standalone function to use numba for efficiency."""
    ham = 0
    for ii in range(0, N - 1):
        for jj in range(ii + 1, N):
            ham = ham + spins[ii] * spins[jj]

    ham = -1 / N * ham
    return ham


def metropolis_step(state_old: State, T: float) -> State:
    """ Performs a MC step in configuration step following the metropolis algorithm."""
    # Choose one spin at random and change it.
    state_new = copy.deepcopy(state_old)  # deepcopy is necessary!!
    state_new.change_spin()

    # Compute difference in hamiltonians
    del_h = state_new.ham - state_old.ham

    # Apply metropolis criterium (in this case it's Boltzmann distribution)
    if del_h < 0:
        state_final = state_new
    else:
        V = np.random.rand()
        if V <= np.exp(-del_h/T):
            state_final = state_new
        else:
            state_final = state_old

    return state_final


def mc_step(state0: State, T: float) -> State:
    """ This function performs a "Monte Carlo step". That is, N "smaller" Metropolis steps."""
    state_current = state0
    for _ in range(0, state0.N):
        state_current = metropolis_step(state_current, T)

    return state_current


def make_plot(stats: Statistics):
    """Plots simulation_results and saves the figure."""
    fig, axs = plt.subplots(2, 2)
    fig.suptitle("Resulting magnitudes as a function of Temperature, N="+str(stats.N), fontweight="bold")

    data = (("Mean magnetization per spin.", "m", stats.mag),
            ("Magnetic susceptibility.", r'$\chi m$', stats.xm),
            ("Total energy.", "E", stats.E),
            ("Heat capacity.", "Cv", stats.cv))

    for ii, ax in enumerate(axs.flat):
        ax.set_xlabel("T")
        ax.set_title(data[ii][0])
        ax.set_ylabel(data[ii][1])
        ax.plot(stats.T_vec, data[ii][2])
        ax.grid()

    plt.tight_layout()
    plt.savefig("./simulation_results/figure_N" + str(stats.N) + ".png", dpi=1200)
    plt.show()


# -----------------------------------------------------------------------------------------------------------
# ----------------- MAIN CODE -------------------------------------------------------------------------------
def main():
    # PARAMS
    N = 10  # Number of particles
    M = 10000  # Number of MonteCarlo steps (for each value of T)
    it_temp = 50  # Number of steps in temperature
    Ti = 0.1
    Tf = 2.5

    # INITIALIZATION
    T_vec = np.linspace(Ti, Tf, it_temp)
    statistics = Statistics(Ti, Tf, it_temp, M, N)

    # ITERATE THROUGH TEMPERATURE.
    for ii in range(it_temp):
        T = float(T_vec[ii])

        # Initialize state and ensemble (one new ensemble for each T)
        ensemble = Ensemble(M, N)
        state = State(N)
        ensemble.add_sample(state)

        # For each T, perform M MC steps
        for jj in range(M-1):
            state = mc_step(state, T)
            ensemble.add_sample(state)

        print('Current T: ', T, ' ', ii+1, '/', it_temp, '\n')
        # os.system('afplay /System/Library/Sounds/Glass.aiff')  # Make a sound once each temperature has finished

        # Compute statistical information
        statistics.average_ensemble(ensemble)

    # COMPUTE FINAL DATA AND SAVE
    statistics.compute_fluctuations()
    statistics.save_file()
    os.system('say "your program has finished"')

    # PLOT
    make_plot(statistics)


if __name__ == "__main__":
    print("Now executing main.\n")
    main()
