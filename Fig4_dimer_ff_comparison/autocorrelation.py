#%%
import pandas as pd
import matplotlib.pyplot as plt
import pymbar.timeseries as ts
import numpy as np
from tqdm import tqdm
import os

# Function to perform analysis on a single CSV file
def analyze_csv(file_path):
    # Load the data
    data = pd.read_csv(file_path)['Interchain SPHF6 Aligned'].values
    print(np.shape(data))
    # look at every 5 frames
    data = data[::5]
    
    # Create a new figure for each file
    plt.figure(figsize=(10, 8))
    
    # plot time series
    plt.subplot(2, 1, 1)
    # on first plot
    plt.plot(data, linewidth=0.1)
    plt.xlim(left=0, right=len(data))
    plt.ylabel('interchain SPHF6 contacts')
    plt.xlabel('frame')
    plt.title(file_path[19:-4])
    
    # Compute the autocorrelation
    C_n = ts.normalized_fluctuation_correlation_function(data)
    
    # Plot the autocorrelation on second plot
    plt.subplot(2, 1, 2)
    plt.plot(C_n)
    plt.axhline(y=0, color='k', linestyle='--')
    plt.xlim(left=0, right=10000)
    plt.ylabel('autocorrelation function')
    plt.xlabel('lag')
    
    # determine when the autocorrelation is less than 0.05
    k_max = np.where(C_n < 0.05)[0]
    if len(k_max) > 0:
        k_max = k_max[0]
    else:
        k_max = len(C_n)
    tau_int = 1 + np.sum(C_n[1:k_max])
    N = len(data)
    N_eff = N / (2 * tau_int)
    
    print(f'Analysis for {file_path}:')
    print('Effective number of samples:', N_eff)
    
    # Save the figure
    plt.tight_layout()
    # plt.savefig(file_path.replace('.csv', '.png'))
    # plt.close()
    plt.show()

# Get list of CSV files in the current working directory
csv_files = [f for f in os.listdir() if f.endswith('.csv')]
# arrange them alphabetically
csv_files.sort()

# Analyze each CSV file
for csv_file in tqdm(csv_files):
    analyze_csv(csv_file)

### storing N_indep values
# C36m:
# a03ws:
# a99SB-disp:

#%%
test_equil_times = np.arange(0,200,10)
num_uncorrelated_samples = []
file_path = 'contacts_HP2_dimer_C36m.csv'
sim_sampling = 2 #ps
my_sampling = 20 #ps
total_sim_time = 300 # ns
integrate_to = 5000 #ps
for equil in tqdm(test_equil_times):
    # Load the data
    colVar = pd.read_csv(file_path)['Interchain SPHF6 Aligned'].values
    skip = round(my_sampling / sim_sampling)
    start = round(equil*1000/sim_sampling)
    colVar = colVar[start::skip]
    C_n = ts.normalized_fluctuation_correlation_function(colVar) # compute the autocorrelation function
    autoC = np.sum(C_n[:round(integrate_to/my_sampling)]) * my_sampling
    num_uncorrelated_samples.append(int((total_sim_time-equil)*1000/(2*autoC)))

plt.plot(test_equil_times,num_uncorrelated_samples,'o')
plt.xlabel('equilibration time (ns)',fontsize=20)
plt.ylabel('number of uncorrelated samples',fontsize=15)

# a99SBdisp: 50ns equil, 578 uncorrelated samples
# a03ws: 50ns equil, 812 uncorrelated samples
# C36m: 0ns equil, 1045 uncorrelated samples (Pritam gave me one with 200ns equil already)