"""This code performs a Monte Carlo simulation of the Ising model considering long-range interactions for a range of T.
Improved version of a code from June 2022, rewritten in January 2024. """
import numpy as np
import os
import numba as nb  # This improves performance, a lot


# -----------------------------------------------------------------------------------------------------------
# ----------------- CLASS DEFINITIONS -----------------------------------------------------------------------
class State:
    """Characterized by the number of particles, their spin (np_array), and the hamiltonian.
    It initializes by default to a random state."""
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

    def change_spin(self, spin_id: int = -1):
        """Changes state of spin, and updates ham. If spin_id is specified (!= -1), it changes said spin.
        Otherwise, by default a random spin is selected."""
        if spin_id == -1:
            spin_id = np.random.randint(0, self.N)
        self.spins[spin_id] *= -1
        self.ham = self.hamiltonian()

    def hamiltonian(self) -> float:
        """Computes the hamiltonian, calling compiled function for efficiency."""
        return hamiltonian_optimized(self.N, self.spins)


class Ensemble:
    """An Ensemble is a collection of M states (each of N particles) used for sampling. Contains a list of lists
    with the spins of each sample, and a list with the hamiltonian of each sample. """
    M: int = 0
    spins: list = []
    ham: list = []

    def add_sample(self, state: State):
        """Adds into the ensemble the desired state, which must be an object of the class State."""
        self.spins.append(state.spins)
        self.ham.append(state.ham)
        self.M += 1


class Statistics:
    """Class to store and compute the statistics of a given Ensemble, for a single value of T, in a list. These include
    the total magnetization M, the mean magnetization per spin m, and the total energy of the state E."""
    N: int
    T: float
    MAG: list
    mag: list
    E: list

    def __init__(self, N: int, T: float):
        self.N = N
        self.T = T
        self.MAG = []
        self.mag = []
        self.E = []

    def add_stats_from_state(self, state: State):
        """Calculates statistics of a single state."""
        MAG_aux = np.abs(np.sum(state.spins))
        ham_aux = state.ham

        self.MAG.append(MAG_aux)
        self.mag.append(1/self.N*MAG_aux)
        self.E.append(ham_aux)

    def save_file(self):
        """Save the final data in a document."""
        filename = "./N"+str(self.N)+"/T"+str(self.T).replace(".", "_")+".txt"
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        save_data = open(filename, "w")
        save_data.write(" ".join([str(element) for element in self.mag]))
        save_data.write("\n")
        save_data.write(" ".join([str(element) for element in self.E]))
        save_data.write("\n")
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


def metropolis_step(state: State, T: float) -> State:
    """Performs a MC step in configuration space following the metropolis algorithm."""
    H0 = state.ham

    # Choose one spin at random, change it and compute new energy:
    spin_id = np.random.randint(0, state.N)
    # state_new = copy.deepcopy(state_old)  # deepcopy is necessary!!
    state.change_spin(spin_id)

    # Compute difference in hamiltonians
    del_h = state.ham - H0

    # Apply metropolis criterium (in this case it's Boltzmann distribution)
    if del_h > 0:
        V = np.random.rand()
        if V >= np.exp(-del_h/T):
            # Change is refused. Change spin again to return to original state.
            state.change_spin(spin_id)

    return state


def mc_step(state: State, T: float) -> State:
    """This function performs a "Monte Carlo step". That is, N "smaller" Metropolis steps."""
    for _ in range(0, state.N):
        state = metropolis_step(state, T)

    return state


# -----------------------------------------------------------------------------------------------------------
# ----------------- MAIN CODE -------------------------------------------------------------------------------
def main():
    # READ PARAMS
    print("Please input parameters\n")
    N = int(input("Number of spins (int): "))
    M = int(input("Number of Monte Carlo steps (int): "))
    it_temp = 50  # Number of steps in temperature
    Ti = 0.1
    Tf = 2.5

    print("\nNOW RUNNING SIMULATION WITH:")
    print(f"Number of spins: {N}")
    print(f"Number of MC steps: {M}")
    print(f"Range of temperatures: Ti={Ti}, Tf={Tf}, steps={it_temp} \n")

    # ITERATE THROUGH TEMPERATURE.
    T_vec = np.linspace(Ti, Tf, it_temp)
    for ii in range(it_temp):

        T = float(T_vec[ii])

        # Initialize state, ensemble and stats (one new for each T)
        statistics = Statistics(N, T)
        ensemble = Ensemble()
        state = State(N)
        ensemble.add_sample(state)

        # For each T, perform M MC steps
        for jj in range(M-1):
            state = mc_step(state, T)
            ensemble.add_sample(state)
            statistics.add_stats_from_state(state)

        statistics.save_file()

        print("\rCurrent T: ", T, " ", ii + 1, "/", it_temp, end="")
        os.system('afplay /System/Library/Sounds/Glass.aiff')  # Make a sound once each temperature has finished

    print("\n\nFinished!")


if __name__ == "__main__":
    main()
