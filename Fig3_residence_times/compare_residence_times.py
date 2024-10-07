#%%
import pandas as pd
import matplotlib.pyplot as plt

# List of residence times files
residence_files = [
    'a03ws_residence_times.csv',
    'a99SBdisp_residence_times.csv',
    'C36m_residence_times.csv'
]

# Dictionary to store residence times data
residence_times_data = {}

# Load residence times data from CSV files
for file in residence_files:
    residence_times_data[file] = pd.read_csv(file).values.flatten().tolist()

# List of colors for each box
colors = ['red', 'blue', 'green']

# Extract data and labels for plotting
data = [residence_times_data[file] for file in residence_files]
labels = ['a03ws', 'a99SBdisp', 'C36m']

# Generate the box plot
plt.figure()

# Create the boxplot with custom colors
boxplot = plt.boxplot(data, patch_artist=True, labels=labels, medianprops=dict(color='black'))

# Color each box
for patch, color in zip(boxplot['boxes'], colors):
    patch.set_facecolor(color)

# Set plot labels and title
plt.ylabel('Residence time (ns)', fontsize=14)
# plt.title('Residence Times for Different Force Fields', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim(bottom=0)
plt.gca().tick_params(axis='x', which='both', bottom=False, top=False)

plt.show()

#%%
import numpy as np

# Define the bin edges spaced every 3 units
bins = np.arange(0, max([max(times) for times in residence_times_data.values()]) + 10, 10)

# Generate the histograms for each dataset
plt.figure(figsize=(10, 6))
for file, color, label in zip(residence_files, colors, labels):
    plt.hist(residence_times_data[file], bins=bins, alpha=0.5, color=color, label=label)

# Set plot labels and title
plt.xlabel('Residence time (ns)', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.title('Histogram of Residence Times for Different Force Fields', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.show()

#%%
import matplotlib.patches as mpatches

# Calculate means and 75th percentiles
means = [np.mean(residence_times_data[file]) for file in residence_files]
percentiles_50 = [np.percentile(residence_times_data[file], 50) for file in residence_files]    
percentiles_90 = [np.percentile(residence_times_data[file], 90) for file in residence_files]

# Labels for the bar graph
labels = ['a03ws', 'a99SBdisp', 'C36m']

# Bar width
bar_width = 0.35

# Positions of the bars on the x-axis
r1 = np.arange(len(labels))/2+ bar_width/2
r2 = [x + bar_width for x in r1]

# List of colors for each force field
colors = ['tab:red', '#039dfc', 'black']

# Create the bar graph
plt.figure(figsize=(5.5, 4.0))

# Bar for 75th percentiles with lower alpha
bars2 = plt.bar(r1, percentiles_90, color=colors, width=bar_width, edgecolor='black', alpha=0.5, label='75th Percentile')

# Bar for means
bars1 = plt.bar(r1, percentiles_50, color=colors, width=bar_width, edgecolor='black', label='Median', hatch='/')

# Apply hatching only to the first two bars
bars1[0].set_hatch('/')
bars1[1].set_hatch('/')
bars1[2].set_hatch('')  # No hatch on the black bar

# Overlay a white-hatched bar on top of the third bar
plt.bar(r1[2], percentiles_50[2], color='none', edgecolor='gray', width=bar_width, hatch='/', label='White Hatch')

# Add labels
plt.xlabel('force field', fontsize=14)
plt.ylabel('residence time (ns)', fontsize=14)
# plt.title('Comparison of Mean and 75th Percentile Residence Times', fontsize=16)
plt.xticks(r1, labels, fontsize=14)
plt.yticks(fontsize=14)
# plt.legend()
# plt.tight_layout()

# redo xticks labels
new_labels = ['a03ws', 'a99SBdisp', 'CHARMM36m']
plt.xticks(r1, new_labels, fontsize=14)

# Define custom patches for the legend
median_patch = mpatches.Patch(facecolor='gray', edgecolor='black', hatch='//', label='median')
percentile_75_patch = mpatches.Patch(facecolor='gray', edgecolor='black', alpha=0.5, label='90th\npercentile')

# Create the custom legend
plt.legend(handles=[median_patch, percentile_75_patch], fontsize=11.5)
# shift up
plt.subplots_adjust(bottom=0.15)
# save figure
plt.savefig('residence_time_comparison.png', dpi=300)

# Show the plot
plt.show()

#%%
