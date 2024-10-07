#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

a03ws_300K = pd.read_csv('fes_a03ws_68CI.csv', index_col=0)
a03ws_325K = pd.read_csv('fes_a03ws_325K_68CI.csv', index_col=0)
a99SBdisp_300K = pd.read_csv('fes_a99SBdisp_68CI.csv', index_col=0)
a99SBdisp_325K = pd.read_csv('fes_a99SBdisp_325K_68CI.csv', index_col=0)
C36m_300K = pd.read_csv('fes_C36m_68CI.csv', index_col=0)
C36m_325K = pd.read_csv('fes_C36m_325K_68CI.csv', index_col=0)

# Put all the dataframes in a list
unprocessed_dfs = [a03ws_300K, a03ws_325K, a99SBdisp_300K, a99SBdisp_325K, C36m_300K, C36m_325K]
dfs = []

# Loop through the dataframes
# Select rows with bincenters between 62.25 and 68.25 (index)
# Then subtract 60.75 from the index
for df in unprocessed_dfs:
    df = df.loc[62.25:68.25]
    df.index -= 60.75
    dfs.append(df)

dG_300K = []
dG_325K = []
# Actually just pull out the last row of each dataframe
for df in dfs[::2]:
    dG_300K.append(df.iloc[-1] * -1)
for df in dfs[1::2]:
    dG_325K.append(df.iloc[-1] * -1)

# Set the width of the bars
bar_width = 0.35

# Extract values and calculate errors for 300K and 325K
dG_300K_vals = [x['f'] for x in dG_300K]
dG_300K_errors = [(x['f_lower'] - x['f'], x['f'] - x['f_upper']) for x in dG_300K]

dG_325K_vals = [x['f'] for x in dG_325K]
dG_325K_errors = [(x['f_lower'] - x['f'], x['f'] - x['f_upper']) for x in dG_325K]

# Calculate error bars
dG_300K_err_lower = [err[0] for err in dG_300K_errors]
dG_300K_err_upper = [err[1] for err in dG_300K_errors]
dG_325K_err_lower = [err[0] for err in dG_325K_errors]
dG_325K_err_upper = [err[1] for err in dG_325K_errors]

# Set positions of bars on X axis
r1 = range(len(dG_300K_vals))
r2 = [x + bar_width for x in r1]

# Colors for each force field
colors_300K = ['tab:red','#039dfc','black']
colors_325K = ['tab:red','#039dfc','black']
hatch_325K = '*'  # Hatch pattern for 325K bars

# Make the plot
plt.figure(figsize=(5.5, 4.0))
ax = plt.gca()

bars1 = ax.bar(r1, dG_300K_vals, color=colors_300K, width=bar_width, edgecolor=(0.2,0.2,0.2),
               yerr=[dG_300K_err_lower, dG_300K_err_upper], capsize=7,ecolor=[0.7,0.7,0.7])
bars2 = ax.bar(r2, dG_325K_vals, color=colors_325K, width=bar_width, edgecolor=(0.2,0.2,0.2),
               yerr=[dG_325K_err_lower, dG_325K_err_upper], capsize=7, hatch=hatch_325K, ecolor=[0.7,0.7,0.7])

# Add labels
ax.set_xlabel('force field    ', fontsize=14)
ax.set_ylabel(r'$ΔG_{dock}$ (kJ/mol)', fontsize=14)
ax.set_xticks([r + bar_width / 2 for r in range(len(dG_300K_vals))])
ax.set_xticklabels(['a03ws', 'a99SBdisp    ', 'CHARMM36m'], fontsize=14)
plt.xticks(fontsize=14), plt.yticks(fontsize=14)

# Define custom patches for the legend
legend_300K = mpatches.Patch(facecolor='white', edgecolor='black', label='300K')
legend_325K = mpatches.Patch(facecolor='white', edgecolor=(0.2,0.2,0.2), hatch='**', label='325K')

# Create the custom legend
plt.legend(handles=[legend_300K, legend_325K], fontsize=14)
# plt.tight_layout()
plt.subplots_adjust(bottom=0.15, left=0.2)
plt.savefig('fes_forcefield_comparison.png', dpi=300)

plt.show()
