#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
# C36m 60ns equilibration time; 1043 samples; integrate 20ns
# a99SBdisp 40ns equilibration time; 179 samples; integrate 20ns
# a03ws 40ns equilibration time; 436 samples; integrate 20ns

# Read and store data for each force field
contacts_list = []
ffs = ['a03ws','a99SBdisp','C36m']
colors = ['tab:red','#039dfc','black']
sampling = [10,10,2] # ps
convergence = [50, 60, 50] # ns
uncorrelated_samples = [179, 1043, 436]

for ff in ffs:
    contacts = pd.read_csv(f'contacts_HP2_dimer_{ff}.csv')['Interchain SPHF6 Aligned']
    equilibrated_contacts = contacts[int(convergence[ffs.index(ff)]*1000/sampling[ffs.index(ff)]):]
    contacts_list.append(equilibrated_contacts)

# Convert list to dictionary for easier access
contacts_dict = dict(zip(ffs, contacts_list))

# Parameters
bootstrap_num = 1000
bin_width = 2.5
bin_range = (0, 60)

# Create bins for histogram
bins = np.arange(bin_range[0], bin_range[1] + bin_width, bin_width)
bin_mids = (bins[1:] + bins[:-1]) / 2

# Bootstrap and calculate probabilities
results = {}
for ff, contacts in contacts_dict.items():
    bootstrapped_histograms = []
    for _ in tqdm(range(bootstrap_num), desc=f'Bootstrapping {ff}'):
        # Sample with replacement
        samples = np.random.choice(contacts, size=uncorrelated_samples[ffs.index(ff)], replace=True)
        # Calculate histogram
        hist, _ = np.histogram(samples, bins=bins, density=True)
        bootstrapped_histograms.append(hist)
    
    # Convert list of histograms to a 2D numpy array
    bootstrapped_histograms = np.array(bootstrapped_histograms)
    # Store results
    results[ff] = bootstrapped_histograms

#%%
# Plotting
plt.figure()
for ff, histograms in results.items():
    # Calculate mean and confidence intervals
    mean_probabilities = np.mean(histograms, axis=0)
    lower_bound = np.percentile(histograms, 5, axis=0)
    upper_bound = np.percentile(histograms, 95, axis=0)
    
    # Convert probabilities to free energies
    kb = 0.0083144621  # kJ/mol/K
    T = 300  # K

    padding = 0.00001
    mean_free_energies = -kb * T * np.log(mean_probabilities + padding)
    lower_free_energies = -kb * T * np.log(lower_bound + padding)
    upper_free_energies = -kb * T * np.log(upper_bound + padding)

    # Normalize so that the minimum free energy is 0 for each set
    min_free_energy = np.min(mean_free_energies)
    mean_free_energies -= min_free_energy
    lower_free_energies -= min_free_energy
    upper_free_energies -= min_free_energy
    
    # Plot
    plt.fill_between(bin_mids/2, lower_free_energies, upper_free_energies, alpha=0.3, color=colors[ffs.index(ff)])
    plt.plot(bin_mids/2, mean_free_energies, label=f'{ff}', color=colors[ffs.index(ff)])

plt.ylim(bottom=0, top=14)
plt.xlim(left=0, right = 22)
# plt.title('Free Energy Landscape', )
plt.xlabel('interchain PHF6 aligned contacts', fontsize=14)
plt.ylabel('free energy (kJ/mol)', fontsize=14)
plt.xticks(np.arange(0, 22, 5), fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14,loc='lower left', ncol=3)
plt.subplots_adjust(bottom=0.15)
plt.savefig('bootstrap_interPHF6_energies.png', dpi=300)

plt.show()

#%%