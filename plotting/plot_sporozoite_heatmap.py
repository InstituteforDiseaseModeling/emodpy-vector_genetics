import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab as pl
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams.update({'font.size': 16})

exp_name = 'sporozoite_reduction'
data_dir = os.path.join(os.path.expanduser('~'), 'Github', 'emodpy-vector_genetics', 'data')
fig_dir = os.path.join(os.path.expanduser('~'), 'Github', 'emodpy-vector_genetics', 'figures')
os.makedirs(fig_dir, exist_ok=True)


def new_colormap():
    bottom = pl.cm.get_cmap('Oranges', 128)
    top = pl.cm.get_cmap('Blues_r', 128)
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                           bottom(np.linspace(0, 1, 128))))
    newcmp = mpl.colors.ListedColormap(newcolors, name='OrangeBlue')

    colormap = newcmp

    return colormap

incidence_file = os.path.join(data_dir, 'sporozoite_reduction_incidence_final.csv')
df_incidence = pd.read_csv(incidence_file)
df_elimination = df_incidence[df_incidence['Baseline'] == 0]
df_incidence = df_incidence[df_incidence['Baseline'] == 1]

larval_capacity = [7.5, 7.75, 8]

fig, axs = plt.subplots(2, 3, figsize=(30, 10))
cbar_ax = fig.add_axes([.92, .53, .01, .35])
cbar_ax2 = fig.add_axes([.92, .11, .01, .35])

for k, l in enumerate(larval_capacity):
    df_inc = df_incidence[df_incidence['Larval_Capacity'] == l]
    df_elim = df_elimination[df_elimination['Larval_Capacity'] == l]
    incidence_reduction = (1 - np.array(df_elim['New Clinical Cases']) / np.array(df_inc['New Clinical Cases'])) * 100
    df_inc['Incidence_reduction'] = incidence_reduction

    # Plot incidence reduction
    index = df_inc['Transmission_To_Human'].unique()
    cols = df_inc['Infected_Progress'].unique()
    B = np.reshape(np.array(df_inc['Incidence_reduction']), (-1, 5))
    labels = np.reshape(np.round(np.array(df_inc['Incidence_reduction'])), (-1, 5))
    df_plot = pd.DataFrame(B, columns=cols, index=index)

    sns.heatmap(B, ax=axs[1, k], cmap=new_colormap(), annot=labels, annot_kws={'fontsize': 16}, fmt='0.2g',
                vmin=0, vmax=100, cbar_ax=None if k else cbar_ax2  , cbar=k == 0,
                cbar_kws={'ticks': [0, 25, 50, 75, 100]})
    axs[1, k].invert_xaxis()

    # Plot elimination
    index = df_elim['Transmission_To_Human'].unique()
    cols = df_elim['Infected_Progress'].unique()
    B = np.reshape(np.array(df_elim['Elimination']*100), (-1, 5))
    labels = np.reshape(np.round(np.array(df_elim['Elimination']*100)), (-1, 5))
    df_plot = pd.DataFrame(B, columns=cols, index=index)
    sns.heatmap(B, ax=axs[0, k], cmap="YlGnBu", annot=labels, annot_kws={'fontsize': 16}, fmt='0.4g',
                vmin=0, vmax=100, cbar_ax=None if k else cbar_ax, cbar=k == 0,
                cbar_kws={'ticks': [0, 25, 50, 75, 100]})
    # axs[0, k].invert_yaxis()
    axs[0, k].invert_xaxis()

    axs[0, 1].set_yticks([])
    axs[1, 1].set_yticks([])
    axs[0, 2].set_yticks([])
    axs[1, 2].set_yticks([])
    axs[0, k].set_xticks([])
    xticklabels = ['%i%%' % (100-i) for i in np.array(cols) * 100]
    yticklabels = ['%i%%' % (100-i) for i in np.array(index) * 100]
    axs[1, k].set_xticklabels(xticklabels)
    axs[1, 0].set_yticklabels(yticklabels)
    axs[0, 0].set_yticklabels(yticklabels)
    axs[1, 0].set_ylabel('Sporozoite reduction')
    axs[0, 0].set_ylabel('Sporozoite reduction')
    axs[1, k].set_xlabel('Delay')

    if k == 0:
        ax2 = axs[0, len(larval_capacity) - 1].twinx()
        ax2.set_ylabel('Elimination\nprobability', color='b', rotation=0, labelpad=50,
                       fontdict={'size': 16})
    else:
        ax2 = axs[1, len(larval_capacity) - 1].twinx()
        ax2.set_ylabel('Incidence\nreduction', color='b', rotation=0, labelpad=50,
                       fontdict={'size': 16})
    ax2.set_yticks([])
    for _, spine in ax2.spines.items():
        spine.set_visible(False)

    axs[0, k].set_title('Annual EIR = %i' % df_inc['Annual EIR'].unique()[0])

fig.tight_layout(rect=[0, 0, .92, 0.96])
plt.savefig(os.path.join(fig_dir, 'heatmap.pdf'))
plt.savefig(os.path.join(fig_dir, 'heatmap.png'))
# plt.show()

# df_incidence