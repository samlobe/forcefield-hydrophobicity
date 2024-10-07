#%%
import pymbar.timeseries as timeseries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

dfs = []
ffs = ['C36m','a03ws','a99SBdisp']
for ff in ffs:
    df = pd.read_csv(f'contacts_HP2_dimer_{ff}.csv')
    dfs.append(df)

# set dfs to a dictionary for easy access
dfs = dict(zip(ffs,dfs))

#%%
data = dfs['a03ws']['Interchain SPHF6 Aligned'] # full dataset starting from 0ns
time_step = 10 #ps between frames that were ouput in the trajectory
my_sampling = 20 #ps between frames that I'll load in this code; must be a multiple of time_step
total_time = (len(data)-1)*time_step/1000 #ns total time of the simulation
# equilibration times to try
equilibration_times = np.arange(40,100,10) # 0-200ns in 20ns increments
# number of uncorrelated samples
uncorrelated_samples = []

for i,equilibration in enumerate(tqdm(equilibration_times)):
    starti = int(equilibration*1000/time_step) # index of the first frame to use
    skip = int(my_sampling/time_step) # number of frames to skip
    colVar = data[starti::skip]
    C_n = timeseries.normalizedFluctuationCorrelationFunction(colVar) # compute the autocorrelation function
    xs = np.arange(len(C_n)) * my_sampling / 1000
    # Decide how much of the autocorrelation function to integrate based on the plot
    plt.plot(xs,C_n,lw=1,label=f'{equilibration} ns equil')
    # plot a horizontal line at y=0
    plt.plot([0,xs[-1]],[0,0],'k--',lw=1)

    # I chose to integrate up to 50ns since that's where C_A is completely at 0 for my data
    integrate_up_to = 20
    integrate_index = np.argmin(np.abs(xs-integrate_up_to))
    autoC = np.sum(C_n[:integrate_index]) * my_sampling # integrate the autocorrelation function; `my_sampling` is width of bins in ps
    uncorrelated_samples.append(int((total_time-equilibration)*1000/(2*autoC))) # number of uncorrelated samples


plt.xlim(left=-1,right=100)
plt.xlabel('time (ns)'); plt.ylabel('$C_A$(t)')
plt.legend(fontsize=5,ncol=2)

# Plot the number of uncorrelated samples as a function of the equilibration time
plt.figure()
plt.plot(equilibration_times,uncorrelated_samples)
plt.xlabel('equilibration time (ns)'); plt.ylabel('number of uncorrelated samples')
plt.ylim(bottom=0)

# print the number of uncorrelated samples for each equilibration time
print('Equilibration time (ns) | Number of uncorrelated samples')
for i,equilibration in enumerate(equilibration_times):
    print(f'{equilibration:24.0f} | {uncorrelated_samples[i]:24.0f}')

#%%
# plot the time series of the data
plt.figure()
plt.plot(np.arange(len(data))*time_step/1000,data,lw=.05)

#%%
# C36m 60ns equilibration time; 1043 samples; integrate 20ns
# a99SBdisp 40ns equilibration time; 179 samples; integrate 20ns
# a99SBdisp 540ns equilibration time; 466 samples; integrate 100ns
# a03ws 50ns equilibration time; 436 samples; integrate 20ns
