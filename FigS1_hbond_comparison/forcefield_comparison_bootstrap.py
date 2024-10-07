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
        for line in islice(f,start,None,skip):
            clean_data.append(line.split()[1])
    clean_data = np.asarray(clean_data)
    clean_data = clean_data.astype(float)
    return clean_data

simulations = ['jR2R3_dimer_a03ws_inter_hbond.xvg','jR2R3_dimer_a99SBdisp_inter_hbond.xvg','jR2R3_dimer_C36m_inter_hbond.xvg',
               'jR2R3_P301L_dimer_a03ws_inter_hbond.xvg','jR2R3_P301L_dimer_a99SBdisp_inter_hbond.xvg','jR2R3_P301L_dimer_C36m_inter_hbond.xvg']
labels = ['jR2R3_a03ws','jR2R3_a99SBdisp','jR2R3_C36m',
          'jR2R3_P301L_a03ws','jR2R3_P301L_a99SBdisp','jR2R3_P301L_C36m']
equilibration = [240,100,40,60,60,0] #ns equilibration time from correlation_time2.py
independent_samples = [494,604,397,288,475,637] # from correlation_time2.py
dts = [10,10,2,10,10,2] # ps between each frame

data = [array_cleaning(simulation,equilibration[i],dts[i],dts[i]) for i,simulation in enumerate(simulations)]
#%% BOOTSTRAP
bootstrap_num = 500
bins = np.arange(26)-0.5 # 0,1,2,3,...24 hbonds
# add a bin for 25+ hbonds
bins = np.append(bins,np.inf)
data_bootstrapped = np.zeros((bootstrap_num,len(bins)-1,len(simulations)))
for i in range(len(simulations)):
    for k in range(bootstrap_num):
        indices = np.random.choice(len(data[i]),independent_samples[i])
        boot_data = data[i][indices]
        data_bootstrapped[k,:,i] = np.histogram(boot_data,bins=bins,density=False)[0]/len(indices)
#%% Sort from greatest to least and pull out the 67% confidence interval
confidence = 0.90
sorted_frac = np.sort(data_bootstrapped,axis=0)
il, im, iu = round((1-confidence)/2*bootstrap_num),round(bootstrap_num/2),round((1-(1-confidence)/2)*bootstrap_num)
low_array, mean_array, up_array = sorted_frac[il,:,:], sorted_frac[im,:,:], sorted_frac[iu,:,:]

#%% output csvs
hbond_nums = np.arange(25).astype(str)
hbond_nums = np.append(hbond_nums,'25+')
confidence_levels = ['5%','50%','95%']
counter = 0
dfs = []
for low_row,mean_row,up_row in zip(low_array.T,mean_array.T,up_array.T):
    df = pd.DataFrame(np.vstack((low_row,mean_row,up_row)).T,
                      columns = confidence_levels, index=hbond_nums)
    dfs.append(df)
    df.to_csv(f'{labels[counter]}_confidence.csv')
    counter = counter + 1

#%% Plot bars with errors for jR2R3 in 3 force fields
width = 0.2
colors = ['tab:red', '#039dfc', 'black']
plt.figure(figsize=(8,4))
for i,label in enumerate(labels[0:3]):
    hbond_confidence = pd.read_csv(f'{label}_confidence.csv',index_col=0)
    bins = df.index
    means = hbond_confidence['50%']
    errors = np.zeros((2,len(bins)))
    errors[0,:] = means - hbond_confidence['5%']
    errors[1,:] = hbond_confidence['95%'] - means
    plt.bar(np.arange(len(bins))+width*i,means*100,yerr=errors*100,align='edge',
            width=width,label=label[6:],color=colors[i])
plt.legend(fontsize=14)
plt.xticks(np.arange(len(bins))+.3,bins,fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(left=-0.5,right=18)
plt.ylim(top=90)
plt.xlabel('intermolecular hydrogen bonds',fontsize=14)
plt.ylabel('probability (%)',fontsize=14)
plt.title('jR2R3',fontsize=20)
plt.subplots_adjust(bottom=0.15)
plt.savefig('jR2R3_interhbonds_3ff.png',dpi=300)
plt.show()
#%% Plot bars with errors for jR2R3 P301L in 3 force fields
width = 0.2
plt.figure(figsize=(8,4))
for i,label in enumerate(labels[3:]):
    # print(i)
    hbond_confidence = pd.read_csv(f'{label}_confidence.csv',index_col=0)
    bins = df.index
    means = hbond_confidence['50%']
    errors = np.zeros((2,len(bins)))
    errors[0,:] = means - hbond_confidence['5%']
    errors[1,:] = hbond_confidence['95%'] - means
    plt.bar(np.arange(len(bins))+width*i,means*100,yerr=errors*100,align='edge',
            width=width,label=label[12:],color=colors[i])
plt.legend(fontsize=14)
plt.xticks(np.arange(len(bins))+.3,bins,fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(left=-0.5,right=18)
plt.ylim(top=90)
plt.xlabel('intermolecular hydrogen bonds',fontsize=14)
plt.ylabel('probability (%)',fontsize=14)
plt.title('jR2R3 P301L',fontsize=20)
plt.subplots_adjust(bottom=0.15)
plt.savefig('jR2R3_P301L_interhbonds_3ff.png',dpi=300)
plt.show()

#%% Plot both in each force field
width = 0.2
for i,label in enumerate(labels[0::3]): # ['jR2R3_a03ws', 'jR2R3_P301L_a03ws']
    print(i)
    hbond_confidence = pd.read_csv(f'{label}_confidence.csv',index_col=0)
    bins = df.index
    means = hbond_confidence['50%']
    errors = np.zeros((2,len(bins)))
    errors[0,:] = means - hbond_confidence['5%']
    errors[1,:] = hbond_confidence['95%'] - means
    plt.bar(np.arange(len(bins))+width*i,means*100,yerr=errors*100,align='edge',
            width=width,label=label[:4])
plt.legend()
plt.xticks(np.arange(len(bins))+.3,bins)
plt.xlim(right=18)
plt.ylim(top=90)
plt.xlabel('Number of Intermolecular Hydrogen Bonds',fontsize=14)
plt.ylabel('% Probability',fontsize=14)
plt.title('in a99SB-disp',fontsize=20)
plt.show()

#%%
width = 0.2
for i,label in enumerate(labels[1::3]): # ['jR2R3_a99SBdisp', 'jR2R3_P301L_a99SBdisp']
    print(i)
    hbond_confidence = pd.read_csv(f'{label}_confidence.csv',index_col=0)
    bins = df.index
    means = hbond_confidence['50%']
    errors = np.zeros((2,len(bins)))
    errors[0,:] = means - hbond_confidence['5%']
    errors[1,:] = hbond_confidence['95%'] - means
    plt.bar(np.arange(len(bins))+width*i,means*100,yerr=errors*100,align='edge',
            width=width,label=label[:4])
plt.legend()
plt.xticks(np.arange(len(bins))+.3,bins)
plt.xlim(right=18)
plt.ylim(top=90)
plt.xlabel('Number of Intermolecular Hydrogen Bonds',fontsize=14)
plt.ylabel('% Probability',fontsize=14)
plt.title('in CHARMM36m',fontsize=20)
plt.show()

#%%
width = 0.2
for i,label in enumerate(labels[2::3]): # ['jR2R3_C36m', 'jR2R3_P301L_C36m']
    print(i)
    hbond_confidence = pd.read_csv(f'{label}_confidence.csv',index_col=0)
    bins = df.index
    means = hbond_confidence['50%']
    errors = np.zeros((2,len(bins)))
    errors[0,:] = means - hbond_confidence['5%']
    errors[1,:] = hbond_confidence['95%'] - means
    plt.bar(np.arange(len(bins))+width*i,means*100,yerr=errors*100,align='edge',
            width=width,label=label[:4])
plt.legend()
plt.xticks(np.arange(len(bins))+.3,bins)
plt.xlim(right=18)
plt.ylim(top=90)
plt.xlabel('Number of Intermolecular Hydrogen Bonds',fontsize=14)
plt.ylabel('% Probability',fontsize=14)
plt.title('in a03ws',fontsize=20)
plt.show()