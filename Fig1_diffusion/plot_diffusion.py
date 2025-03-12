#%%
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

data = pd.read_csv('norm_diffusion.csv')

# Extracting the data
amino_acids = data['Amino Acid']
a03wsFF_data = data['a03wsFF_data']
a03wsFF_error = data['a03wsFF_error']
adispFF_data = data['adispFF_data']
adispFF_error = data['adispFF_error']
charmm36mFF_data = data['charmm36mFF_data']
charmm36mFF_error = data['charmm36mFF_error']

# Use 1-letter code for amino acids
three_to_one = {'Ala':'A', 'Arg':'R', 'Asn':'N', 'Asp':'D', 'Cys':'C', 'Gln':'Q', 'Glu':'E', 'Gly':'G', 'His':'H', 'Ile':'I', 'Leu':'L', 'Lys':'K', 'Met':'M', 'Phe':'F', 'Pro':'P', 'Ser':'S', 'Thr':'T', 'Trp':'W', 'Tyr':'Y', 'Val':'V'}
amino_acids = [three_to_one[aa] for aa in amino_acids]

plt.figure(figsize=(10, 4))
# Creating the scatter plot
colors = ['tab:red', '#039dfc', 'black']
plt.errorbar(amino_acids, a03wsFF_data, yerr=a03wsFF_error, fmt='o', label='a03ws', color=colors[0], capsize=5)
plt.errorbar(amino_acids, adispFF_data, yerr=adispFF_error, fmt='o', label='a99SBdisp', color=colors[1], capsize=5)
plt.errorbar(amino_acids, charmm36mFF_data, yerr=charmm36mFF_error, fmt='o', label='CHARMM36m', color=colors[2], capsize=5)

plt.xlabel('amino acid',fontsize=15)
plt.ylim(0.6,1)
plt.ylabel(r'normalized diffusion $(D/D_{bulk})$',fontsize=16)
# plt.title('Normalized Diffusion Data by Force Field',fontsize=12)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16, loc='center', bbox_to_anchor=(0.2, 0.62))

plt.yticks([0.6, 0.7, 0.8, 0.9, 1.0], fontsize=16)
plt.subplots_adjust(bottom=0.15)

# plt.tight_layout()
# Save the plot
plt.savefig('norm_diffusion.png', dpi=300)

# Display the plot
plt.show()

#%%
import numpy as np
import matplotlib.pyplot as plt

def custom_transform(y, cutoff=0.90, stretch=4):
    # Convert y to a NumPy array for element-wise operations
    y = np.array(y)
    # Manually stretch the range from the cutoff to 1.0
    stretched_y = np.where(y > cutoff, stretch * (y - cutoff) + cutoff, y)
    return stretched_y

def inverse_transform(y, cutoff=0.90, stretch=4):
    # Reverse the custom transformation
    original_y = np.where(y > cutoff, (y - cutoff) / stretch + cutoff, y)
    return original_y

cutoff=0.9
stretch = 3

plt.figure(figsize=(10, 4))
# Creating the scatter plot
colors = ['tab:red', '#039dfc', 'black']
plt.errorbar(amino_acids, custom_transform(a03wsFF_data, cutoff=cutoff, stretch=stretch), yerr=a03wsFF_error, fmt='o', label='a03ws', color=colors[0], capsize=5)
plt.errorbar(amino_acids, custom_transform(adispFF_data, cutoff=cutoff, stretch=stretch), yerr=adispFF_error, fmt='o', label='a99SBdisp', color=colors[1], capsize=5)
plt.errorbar(amino_acids, custom_transform(charmm36mFF_data, cutoff=cutoff, stretch=stretch), yerr=charmm36mFF_error, fmt='o', label='CHARMM36m', color=colors[2], capsize=5)

plt.xlabel('amino acid', fontsize=15)
plt.ylabel(r'normalized diffusion $(D/D_{bulk})$', fontsize=16)

# Add shaded region for stretched area

plt.axhspan(custom_transform(cutoff), custom_transform(1.0,cutoff=cutoff,stretch=stretch), color=(0.8,0.8,0.8), alpha=0.2)
plt.ylim(0.6, 1+0.1*(stretch-1))

# Set custom y-ticks with inverse transformation
original_ticks = [0.6, 0.7, 0.8, 0.9, 0.95, 1.0]  # Include extra ticks in the stretched range
plt.yticks(custom_transform(original_ticks, cutoff=cutoff, stretch=stretch), [f'{tick:.2f}' for tick in original_ticks], fontsize=16)
plt.xticks(fontsize=16)
plt.axhline(y=0.9, color=[0.9,0.9,0.9], linestyle='--', lw=1, alpha=1)
plt.axhline(y=0.8, color=[0.9,0.9,0.9], linestyle='--', lw=1, alpha=1)
plt.axhline(y=0.7, color=[0.9,0.9,0.9], linestyle='--', lw=1, alpha=1)

# add annotated labels
plt.text(-0.5, 1.13, 'CHARMM36m', color='black', fontsize=16)
plt.text(-0.5, 0.93, 'a03ws', color='tab:red', fontsize=16)
plt.text(-0.5, 0.75, 'a99SBdisp', color='#039dfc', fontsize=16)

# plt.legend(fontsize=16, loc='center', bbox_to_anchor=(0.2, 0.42))
plt.subplots_adjust(bottom=0.15)
# save the plot
plt.savefig('stretched_norm_diffusion.png', dpi=300)
plt.show()

