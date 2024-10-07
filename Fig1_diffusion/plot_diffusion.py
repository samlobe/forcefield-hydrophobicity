#%%
import pandas as pd
import matplotlib.pyplot as plt

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
# shift up slightly for better visibility
plt.subplots_adjust(bottom=0.15)
# plt.tight_layout()
# Save the plot
plt.savefig('norm_diffusion.png', dpi=300)

# Display the plot
plt.show()

