#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ffs = ['a03ws','a99SBdisp','C36m']
dfs = []
for ff in ffs:
    df = pd.read_csv(f'obs_box_3body_{ff}.csv',index_col=0)
    dfs.append(df)

# put into dictionary
data = {}
for i,ff in enumerate(ffs):
    data[ff] = dfs[i]

#%%
# measure icosahedral fraction
processed_data = {}
for ff in ffs:
    df = data[ff]
    # definition: 60-65 degrees
    icosa = df.loc[:,'62.5']*5 # 5 degree bin width
    # divide each row by the bulk value
    rel_icosa = icosa/icosa['bulk']
    # subtract each row by the bulk values
    diff = df.sub(df.loc['bulk'],axis=1)
    icosa_diff = diff.loc[:,'62.5']*5 # 5 degree bin width
    # save new dataframe with tetra and tetra_frac
    processed_data[ff] = pd.DataFrame({'icosa':icosa,'rel_icosa':rel_icosa,'icosa_diff':icosa_diff})

#%%
colors = ['tab:red','#039dfc','black']
fig, ax = plt.subplots(figsize=(2.75,4))
for i,ff in enumerate(ffs):
    df = processed_data[ff]
    selected_data = [df.loc['y=20-24', 'rel_icosa'], df.loc['y=5-10', 'rel_icosa']]
    percent_inc = 100*(np.array(selected_data)-1)
    print(ff,percent_inc)
    ax.plot([0,1],selected_data,marker='s',label=ff,color=colors[i])
# set y-ticks as 'y=5-10 and y=20-24'
ax.set_xticks([0,1])
ax.set_xticklabels(['undocked\npeptide\n(14<y<18Å)','docked\npeptide\n(y < 4Å)'],fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel("relative water icosahedrality\nin dashed box (global)",fontsize=12)
# plt.legend(loc='upper center' ,fontsize=12)
plt.text(0.5,0.979,'a03ws',fontsize=12,ha='center',color=colors[0],rotation=16)
plt.text(0.5,0.963,'a99SBdisp',fontsize=12,ha='center',color=colors[1],rotation=35)
plt.text(0.5,0.991,'CHARMM36m',fontsize=12,ha='center',color=colors[2],rotation=9)
ax.yaxis.set_major_formatter(plt.matplotlib.ticker.StrMethodFormatter('{x:,.2f}'))
plt.ylim(bottom=0.95,top=1.00)
plt.tight_layout()
plt.savefig('obs_box_3body_rel_icosa.png',dpi=300)
