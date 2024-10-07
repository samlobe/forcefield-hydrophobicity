#%%
import pandas as pd
import matplotlib.pyplot as plt

data_300K = pd.read_csv('fes_a99SBdisp_68CI.csv',index_col=0)
data_325K = pd.read_csv('fes_a99SBdisp_325K_68CI.csv',index_col=0)

# get just the rows with between bincenters 62.25 and 68.25 (index)
data_300K = data_300K.loc[62.25:68.25]
data_325K = data_325K.loc[62.25:68.25]

plt.figure(figsize=(5.5, 4.0))

# subtract 60.75 from the index
data_300K.index -= 60.75
data_325K.index -= 60.75

# Compute the error margins
yerr_lower = data_300K['f'] - data_300K['f_lower']
yerr_upper = data_300K['f_upper'] - data_300K['f']

# Plot the profile, shading the 90% confidence interval
plt.plot(data_300K.index, data_300K['f'], 'o-', label='a99SBdisp: T=300K',
         color='#039dfc', markersize=8,markeredgecolor='black', markeredgewidth=1)
plt.fill_between(data_300K.index, data_300K['f_lower'], data_300K['f_upper'],
                 color='#039dfc', alpha=0.5)

# repeat for 325K data
yerr_lower = data_325K['f'] - data_325K['f_lower']
yerr_upper = data_325K['f_upper'] - data_325K['f']

plt.plot(data_325K.index, data_325K['f'], '*-', label='a99SBdisp: T=325K',
         color='#9828a6', markersize=10, markeredgecolor='black', markeredgewidth=1)
plt.fill_between(data_325K.index, data_325K['f_lower'], data_325K['f_upper'],
                 color='#9828a6', alpha=0.5)

plt.xlabel('center-of-mass y-coordinate (Å)', fontsize=14)
plt.ylabel('free energy (kJ/mol)', fontsize=14)
# plt.title('a99SB-disp: NaCl 100mM', fontsize=16)
plt.xticks(fontsize=14), plt.yticks(fontsize=14)
plt.ylim(0,2)
plt.xlim(1.5,8.5)
# only include x-ticks from 2 to 7
plt.xticks(range(2,9))
# include tick marks on the right (without numbers)
# plt.tick_params(right=True)

plt.legend(loc='upper left', fontsize=14)
plt.tight_layout()
plt.savefig('fes_a99SBdisp_300-325K.png', dpi=300)

plt.show()