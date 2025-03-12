#%%
import pymbar.timeseries as timeseries
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import islice

#%%
ffs = ['a99SBdisp','C36m','a03ws']
labels = ['a99SBdisp','C36m','a03ws']
files = {'a99SBdisp':'com_distances_a99SBdisp.csv',
         'C36m':'com_distances_C36m.csv',
         'a03ws':'com_distances_a03ws.csv'}
sampling_freqs = {'a99SBdisp':10,'C36m':10,'a03ws':10}
# load the number of intermolecular hydrogen bond per frame (created with  `gmx hbond ...` with mainchain+H selected for each monomer)
data = []
for ff in ffs:
    data.append(np.loadtxt(f'com_distances_{ff}.csv',skiprows=1))

#%%
# df = pd.DataFrame(np.vstack((data,list(sampling_freqs.values()))),
#                   index=['com_distance','sampling_freq'],
#                   columns=labels)
df = {}
for i,ff in enumerate(ffs):
    df[ff] = {'com_distance':data[i],'sampling_freq':sampling_freqs[ff]}


#%%
# this will pull out the 2nd column of your textfile. There can't be any comments or it wont work.
def array_cleaning(file_name, start_time, simulation_sampling_freq, my_sampling_freq):
    # start_time in ns
    # simulation_sampling_freq = the time step of the data in the file (ps)
    # my_sampling_freq = the time step you want (ps); maybe 2x the autocorrelation time
    skip = round(my_sampling_freq / simulation_sampling_freq)
    start = round(start_time*1000/simulation_sampling_freq)
    clean_data = []
    with open(file_name, 'r') as f:
        for line in islice(f,start+1,None,skip):
            # print(line)
            clean_data.append(line.split()[0])
    clean_data = np.asarray(clean_data)
    clean_data = clean_data.astype(float)
    return clean_data

calc_autoc_in_block = True # to print the autocorrelation time in the following block of your simulation
block_start = 100 #ns
block_end = 500 #ns
plot_autoc_over_time = True # to plot how the autocorrelation time changes throughout the simulation

### LOAD DATA ###
which_ff = 'C36m'
file_name = files[which_ff]
sim_sampling =  sampling_freqs[which_ff] #ps frequency of structures in my simulation data
my_sampling = 20 #ps frequency of structures I'll load in this code
colVar = array_cleaning(files[which_ff],0,sim_sampling,my_sampling)
print(f'Looking at the variable in {files[which_ff]}:\n')

## CALCULATE AUTOCORRELATION TIME OVER THE FULL SIMULATION (autoc_time = 1/2 statistical inefficiency)
full_autoc_time = my_sampling * timeseries.integrated_autocorrelation_time(colVar)
print(f'The autocorrelation time over the FULL simulation is {full_autoc_time:.0f} ps')

## CALCULATE AUTOCORRELATION TIME FROM <block_start> ns to <block end> ns
if calc_autoc_in_block == True:
    xx = block_start * 1000 / my_sampling
    yy = block_end   * 1000 / my_sampling
    autoc_from_xx_to_yy = my_sampling * timeseries.integrated_autocorrelation_time(colVar[int(xx):int(yy)])
    print(f'The autocorrelation time from {block_start}-{block_end}ns is {autoc_from_xx_to_yy:.0f} ps')

## SEE HOW AUTOCORRELATION TIME CHANGES THROUGHOUT THE SIMULATION
# Each data point represents the autocorrelation time from that time point going forward 100ns
if plot_autoc_over_time == True:
    fig = plt.figure(figsize=(10, 8)) # create figure object with our desired size
    
    end_time = len(colVar) * my_sampling / 1000 # end time of the simulation in ns
    spacing = 10 # ns between the data I want to plot
    times = np.arange(0,end_time - 100 - spacing, spacing) # 0, 10, 20 ... 500 if end_time is 600 and spacing = 10.
    
    autoc_100ns_block = []
    for time in times:
        begin = time * 1000 / my_sampling # index of data to begin
        end = begin + (100 * 1000) / my_sampling # index of data to end (100ns after begin)
        autoc_100ns_block.append(my_sampling * timeseries.integrated_autocorrelation_time(colVar[int(begin):int(end)])) # calculates autocorrelation time for that 100ns block
    plt.plot(times,autoc_100ns_block,'o',label = file_name)
        
    plt.xlabel('time (ns)',fontsize=20)
    plt.ylabel('autocorrelation time (ps)',fontsize=20)
    plt.title('How does autocorrelation time (looking forward 100ns) \n change throughout the simulation?',fontsize=20)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.legend(prop={'size': 15})
    
    plt.show()
#%% Now let's find the autocorrelation times
which_data = 'a03ws' # swap in each force field here
data = df[which_data]['com_distance'] # full dataset starting from 0ns
time_step = df[which_data]['sampling_freq'] #ps between frames that were ouput in the trajectory
my_sampling = 20 #ps between frames that I'll load in this code; must be a multiple of time_step
total_time = (len(data)-1)*time_step/1000 #ns total time of the simulation
# equilibration times to try
equilibration_times = np.arange(0,200,10) # 0-200ns in 20ns increments
# number of uncorrelated samples
uncorrelated_samples = []

for i,equilibration in enumerate(tqdm(equilibration_times)):
    starti = int(equilibration*1000/time_step) # index of the first frame to use
    skip = int(my_sampling/time_step) # number of frames to skip
    colVar = data[starti::skip]
    C_n = timeseries.normalized_fluctuation_correlation_function(colVar) # compute the autocorrelation function
    xs = np.arange(len(C_n)) * my_sampling / 1000
    # Decide how much of the autocorrelation function to integrate based on the plot
    plt.plot(xs,C_n,lw=1,label=f'{equilibration} ns equil')
    # plot a horizontal line at y=0
    plt.plot([0,xs[-1]],[0,0],'k--',lw=1)

    # I chose to integrate up to 50ns since that's where C_A is completely at 0
    integrate_up_to = 10
    integrate_index = np.argmin(np.abs(xs-integrate_up_to))
    autoC = np.sum(C_n[:integrate_index]) * my_sampling # integrate the autocorrelation function; `my_sampling` is the width of each bin in ps
    uncorrelated_samples.append(int((total_time-equilibration)*1000/(2*autoC))) # number of uncorrelated samples

plt.xlim(left=-1,right=50); plt.ylim(top=1)
plt.xlabel('time (ns)'); plt.ylabel('$C_A$(t)')
plt.legend(fontsize=5,ncol=2)

# Plot the number of uncorrelated samples as a function of the equilibration time
plt.figure()
plt.plot(equilibration_times,uncorrelated_samples)
plt.xlabel('equilibration time (ns)'); plt.ylabel('number of uncorrelated samples')

#%% Manually pick the start time from the plot above
asserted_equilibration = 40 # ns
uncorrelated_samples = uncorrelated_samples[np.argmin(np.abs(equilibration_times-asserted_equilibration))] #141 uncorrelated samples
print(f'{uncorrelated_samples} uncorrelated samples for {which_data}')

### NOTES
# SLS2_a99SB-disp:  40ns equilibration, 667 uncorrelated samples
# SLS2_CHARMM36m:   60ns equilibration, 1234 uncorrelated samples
# SLS2_a03ws:       40ns equilibration, 495 uncorrelated samples


# %%
