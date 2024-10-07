#%%
import pandas as pd
import matplotlib.pyplot as plt

# Load the new data from the uploaded file
file_path = 'stacked_free_energies.csv'
data = pd.read_csv(file_path)

# Convert the specified order to uppercase
amino_acids_order = ['Leu', 'Phe', 'Ile', 'Val', 'Trp', 'Tyr', 'Met', 'Ala', 'Pro', 'Cys', 'Arg', 'His', 'Gly', 'Thr', 'Asn', 'Ser', 'Lys', 'Gln', 'Asp', 'Glu']
amino_acids_order_upper = [aa.upper() for aa in amino_acids_order]

# Pivot the data for easier plotting
pivot_data = data.pivot(index='AA', columns='FF', values='free_energy')

# Reindex the pivot data to match the specified order
pivot_data = pivot_data.reindex(amino_acids_order_upper)

# rename the AA to the one-letter code
# Convert to single letter code
three_to_one = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
                'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I',
                'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
                'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
amino_acids_order_upper = [three_to_one[aa] for aa in amino_acids_order_upper]
pivot_data.index = amino_acids_order_upper

# reorder columns to 'a03ws', 'a99SBdisp', 'C36m'
pivot_data = pivot_data[['a03ws', 'a99SBdisp', 'C36m']]

#%%
# Create a bar plot for the data with adjusted bar width and spacing between groups
colors = ['#B63841','#006B9F','#FEDD4D']
colors = ['#D55E00','#3070AD','#EEE461']
colors = ['tab:red', 'tab:blue', 'gold']
colors = ['tab:red', 'tab:blue', 'black']
colors = ['tab:red','#039dfc','black']
# plt.figure(figsize=(8, 4))
fig,ax = plt.subplots(figsize=(10, 4))
pivot_data.plot(ax=ax, kind='bar', color=colors, width=0.7, position=0.5, edgecolor='black')

plt.xlabel('amino acid', fontsize=16)
plt.xticks(rotation=0, fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel(r'$F_{dewet}$ (kJ/mol per water)', fontsize=16)
plt.legend(['a03ws', 'a99SBdisp', 'CHARMM36m'], loc='upper left', fontsize=16)
plt.ylim(2,6)
# shift up slightly for better visibility
plt.subplots_adjust(bottom=0.15)
# plt.tight_layout()
# Save the plot
plt.savefig('Fdewet.png', dpi=300)

# Display the plot
plt.show()


#%%
# # Create a bar plot for the data
# pivot_data.plot(kind='bar', color=['red', 'blue', 'green'], width=0.8)

# plt.xlabel('Amino Acid', fontsize=12)
# plt.xticks(rotation=0, fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylabel(r'$F_{dewet}$ (kJ/mol per water)', fontsize=12)
# plt.legend(['a03ws', 'a99SBdisp', 'CHARMM36m'], loc='upper left', fontsize=12)
# plt.ylim(2,6)
# # plt.tight_layout()
# # Save the plot
# plt.savefig('Fdewet.png', dpi=300)

# # Display the plot
# plt.show()
