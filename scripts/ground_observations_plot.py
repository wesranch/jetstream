import matplotlib.pyplot as plt
import pandas as pd
import glob

"""

author: @wesranch

Boxplot of CAFI data across the plot network 1994-2017
Some years have sparse observations
    
"""

path = "H:/MS-Research/alaska/data/CooperativeAlaska/biomass*"
files = glob.glob(path)

# clean names for plot    
dfs = list(map(pd.read_csv, files))
df = pd.concat(dfs, axis=0)
df = df[df['Species']!= 'black cottonwood']


#boxplot
species_list = ['resin birch','black spruce', 'white spruce', 'quaking aspen']
carbon_df = [df[df['Species'] == species]['Bio.g.m2'] * 0.47 for species in species_list]

labels = ['Alaskan birch','Black spruce', 'White spruce', 'Quaking aspen']
colors = ['#67161CFF', '#3F6148FF', '#A4804CFF', '#4B5F80FF']
species_counts = [df[df['Species'] == species].shape[0] for species in species_list]

fig, ax = plt.subplots()
ax.set_ylabel('Aboveground C (g/m2)')

bplot = ax.boxplot(carbon_df,
                   patch_artist=True,
                   tick_labels=labels,
                   medianprops=dict(color='black'))
for i, count in enumerate(species_counts):
    ax.text(i + 1, max(carbon_df[i]) + 10, f'# obs. {count}', 
            horizontalalignment='center', 
            fontsize=10, color='black', verticalalignment='bottom')
    
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)

plt.savefig('H:/MS-Research/alaska/SPRING_GRAPHS/CarbonDistributionGround.png', transparent = False, dpi = 300)
plt.show()
