#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import islice

def array_cleaning(file_name, start_time, simulation_sampling_freq, my_sampling_freq):
    # start_time in ns
    # simulation_sampling_freq = the time step of the data in the file (ps)
    # my_sampling_freq = the time step you want (ps); maybe 2x the autocorrelation time
    skip = round(my_sampling_freq / simulation_sampling_freq)
    start = round(start_time*1000/simulation_sampling_freq)
    clean_data = []
    with open(file_name, 'r') as f:
        for line in islice(f,start+1,None,skip): # +1 ignores the header
            clean_data.append(line)
    clean_data = np.asarray(clean_data)
    clean_data = clean_data.astype(float)
    return clean_data

simulations = ['a99SBdisp','C36m','a03ws']
labels = ['a99SBdisp','C36m','a03ws']
equilibration = [40,60,40] #ns equilibration time from correlation_time2.py
independent_samples = [667,1234,495] # from correlation_time2.py
dts = [10,10,10] # ps between each frame

data = [array_cleaning(f'com_distances_{simulation}.csv',equilibration[i],dts[i],dts[i]) for i,simulation in enumerate(simulations)]
#%% BOOTSTRAP
bootstrap_num = 500
bins = np.linspace(0,5,11)
data_bootstrapped = np.zeros((bootstrap_num,len(bins)-1,len(simulations)))
for i in range(len(simulations)):
    for k in range(bootstrap_num):
        indices = np.random.choice(len(data[i]),independent_samples[i])
        boot_data = data[i][indices]
        data_bootstrapped[k,:,i] = np.histogram(boot_data,bins=bins,density=False)[0]/len(indices)
#%% Sort from greatest to least and pull out the 90% confidence interval
confidence = 0.90
sorted_frac = np.sort(data_bootstrapped,axis=0)
il, im, iu = round((1-confidence)/2*bootstrap_num),round(bootstrap_num/2),round((1-(1-confidence)/2)*bootstrap_num)
low_array, mean_array, up_array = sorted_frac[il,:,:], sorted_frac[im,:,:], sorted_frac[iu,:,:]

#%% output csvs
bin_mids = (bins[1:]+bins[:-1])/2 # middle of bins in nm
bin_mid_labels = [f'{bin_mid} nm' for bin_mid in bin_mids]
confidence_levels = ['5%','50%','95%']
counter = 0
dfs = []
for low_row,mean_row,up_row in zip(low_array.T,mean_array.T,up_array.T):
    df = pd.DataFrame(np.vstack((low_row,mean_row,up_row)).T,
                      columns = confidence_levels, index=bin_mid_labels)
    dfs.append(df)
    df.to_csv(f'{labels[counter]}_confidence.csv')
    counter = counter + 1

# plotting the edge of the first and last bin
bin_mids = np.append(0,bin_mids)
bin_mids = np.append(bin_mids,5)


#%% for final figure
fig,ax=plt.subplots(figsize=(6,4))
colors = ['tab:red','#039dfc','black']
for i,label in enumerate(['a03ws','a99SBdisp','C36m']):
    print(label)
    com_confidence = pd.read_csv(f'{label}_confidence.csv',index_col=0)
    low = np.append(0,com_confidence['5%'].values)
    low = np.append(low,0)
    mean = np.append(0,com_confidence['50%'].values)
    mean = np.append(mean,0)
    high = np.append(0,com_confidence['95%'].values)
    high = np.append(high,0)
    #bins = df.index
    # Plot line plot of means
    if label == 'C36m': label ='CHARMM36m'
    plt.plot(bin_mids, mean, label=label,color=colors[i])
    # Plot shaded error bounds
    plt.fill_between(bin_mids, low, high, alpha=0.3,color=colors[i])

plt.legend(fontsize=14)
plt.ylim([0,0.55])
plt.xlim([0,5])
plt.xlabel('dimer center-of-mass distance (nm)',fontsize=14)
plt.ylabel('probability density',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.subplots_adjust(bottom=0.15)
plt.savefig('dimer_com.png',dpi=300)
plt.show()