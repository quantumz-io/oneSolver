
__copyright__ = "hyperQ – Ewa Hendzel"
__license__ = "MIT"
__version__ = "0.1.0"
__status__ = "Dev"


import matplotlib.pylab as plt
import matplotlib as mpl
import dimod
import pandas as pd
import seaborn as sns

plt.style.use('fivethirtyeight')

SIZES = [128, 512, 2048]

for size in SIZES:
    result_df = pd.read_csv(f"./results/{size}/001.csv")
    fig, axes = plt.subplots(1, 3, sharey=False, figsize=(20, 8.), gridspec_kw={"wspace": 0.4})
    
    best_df = result_df.query("num_tries == 460")
    half_tries_df = result_df.query("num_tries == 260")
    
    with open(f"results/{size}/groundstates_TN.txt") as ground_state_file:
        ground_line = list(ground_state_file)[0].strip()
    
    with open(f"chimera_droplets/{size}power/001.txt") as bqm_file:
        bqm = dimod.BQM.from_coo(bqm_file, vartype="SPIN")
        
    bqm.change_vartype("BINARY", inplace=True)
    bqm.offset = 0 # Because offset is now nonzero due to conversion Ising -> Qubo
    
    sample = {
        i: int(bit) for i, bit in enumerate(ground_line[ground_line.find(":")+2:].split(" ")[1:], start=1)
    }
        
    ground_energy = bqm.energy(sample)
    
    for num_tries, ax in zip([160, 310, 460], axes):
        df = result_df.query(f"num_tries == {num_tries}")
        linear_df = df.query("schedule == 'linear'")
        geometric_df = df.query("schedule == 'geometric'")

        linear_energies = linear_df.groupby("num_iter").energy.min() - ground_energy
        geometric_energies = geometric_df.groupby("num_iter").energy.min() - ground_energy

        linear_energies.plot(ax=ax, label="linear", c="red")
        geometric_energies.plot(ax=ax, label="geometric", c="blue")

        ax.set_ylabel("Distance from ground", fontweight ='bold', labelpad=30.)
        ax.set_xlabel("Annealing time", fontweight ='bold', labelpad=5.)
        ax.set_title(f"Trajectories: {num_tries}")
  
    suptitle = fig.suptitle(f"HYP-24 & HYP 23: benchmarks of one-solver-anneal (Chimera of size {size})")
    suptitle.set_position((0.5, 1.05))

    plt.legend()
    plt.savefig(f"HYP-24_chimera_{size}.png")

