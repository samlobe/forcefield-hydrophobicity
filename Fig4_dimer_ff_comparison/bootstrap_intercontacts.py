#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from autocorrelation.py
# a99SBdisp: 50ns equil, 578 uncorrelated samples
# a03ws: 50ns equil, 812 uncorrelated samples
# C36m: 0ns equil, 1045 uncorrelated samples (Pritam gave me one with 200ns equil already)
equilibration_dict = {'a99SBdisp': 50, 'a03ws': 50, 'C36m': 50} # ns
uncorrelated_samples_dict = {'a99SBdisp': 578, 'a03ws': 812, 'C36m': 841}
sampling_dict = {'a99SBdisp': 10, 'a03ws': 10, 'C36m': 10} # ps
ffs = ['a03ws','a99SBdisp','C36m']
confidence_intervals = {}

for ff in ffs:
    file_path = f'contacts_HP2_dimer_{ff}.csv'
    data = pd.read_csv(file_path)['Interchain SPHF6 Aligned'].values / 2
    start = round(equilibration_dict[ff]*1000/sampling_dict[ff])
    data = data[start:]
    # histogram the data
    plt.hist(data, bins=100, density=True)
    plt.title(ff)
    plt.xlabel('number of intercontacts')
    plt.ylabel('probability density')
    plt.show()

    bootstrap = 1000
    n_samples = uncorrelated_samples_dict[ff]
    intercontacts_frac = []
    cutoff = 12 # looking at fraction above this num of intercontacts

    for i in range(bootstrap):
        sample = np.random.choice(data, n_samples)
        intercontacts_frac.append(np.sum(sample > cutoff)/len(sample))

    intercontacts_frac = np.array(intercontacts_frac)
    # get mean, 5th, and 95th percentiles
    mean = np.mean(intercontacts_frac)
    low = np.percentile(intercontacts_frac, 5)
    high = np.percentile(intercontacts_frac, 95)

    # store in confidence_intervals
    confidence_intervals[ff] = [mean, low, high]

#%% Plot bars with error bars
fig, ax = plt.subplots(figsize=(5,4))

means = np.array([confidence_intervals[ff][0] for ff in ffs])*100
low = np.array([confidence_intervals[ff][1] for ff in ffs])*100
high = np.array([confidence_intervals[ff][2] for ff in ffs])*100

# put means, low, high in a csv
df = pd.DataFrame({'FF':ffs,'mean':means,'low':low,'high':high})
df.to_csv('intercontacts_bootstrap.csv',index=False)

colors = ['tab:red', '#039dfc', 'black']
x = np.arange(len(ffs))
ax.bar(x, means, yerr=[means-low, high-means], align='center', alpha=1, ecolor='gray', color=colors,capsize=5)
ax.set_xticks(x)
fontsize=14
ax.set_xticklabels(ffs,fontsize=fontsize)
ax.set_xticklabels(['a03ws','a99SBdisp','CHARMM36m'],fontsize=fontsize)
plt.xlabel('force field',fontsize=fontsize)
plt.yticks(fontsize=fontsize)
ax.set_ylabel('% probability of >12 interchain\naligned PHF6 contacts',fontsize=fontsize)
# save fig
plt.tight_layout()
plt.savefig('intercontacts_bootstrap.png',dpi=300)
plt.show()



# %%
