import pandas as pd
import numpy as np
import os
import math
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

prev_file = os.path.join(data_dir, 'sporozoite_reduction_prevalence_final.csv')
allele_file = os.path.join(data_dir, 'sporozoite_reduction_allele_freqs.csv')

def plot_prevalence():
    df = pd.read_csv(prev_file)
    df = df[df['Transmission_To_Human'] > 0.1]
    df_baseline = df[df['Baseline'] == 1]
    df_baseline = df_baseline.groupby(['Larval_Capacity', 'Time'])[['True Prevalence',
                                                                    'True Prevalence_std']].apply(np.mean).reset_index()
    df = df[df['Baseline'] == 0]

    lhabitat = list(df['Larval_Capacity'].unique())
    tran_human = list(df['Transmission_To_Human'].unique())
    EIR = [15, 30, 60]
    fig_labels = dict(zip(lhabitat, EIR))

    column = 'True Prevalence'
    majorLocator = MultipleLocator(365)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(365 / 12.0)

    cmap = mpl.cm.Blues
    # extract all colors from the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]

    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'Custom cmap', cmaplist[20:], len(tran_human))

    # define the bins and normalize
    bounds = np.linspace(0, len(tran_human), len(tran_human)+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    for i, (g, df_g) in enumerate(df.groupby('Larval_Capacity')):
        fig, axs = plt.subplots(2, 3, figsize=(15, 8))
        df_b = df_baseline[df_baseline['Larval_Capacity']==g]
        for j, (delay, df_delay) in enumerate(df_g.groupby('Infected_Progress')):
            x = math.floor(j / 3)
            y = j % 3

            axs[x, y].plot(df_b['Time'], df_b[column], color='gray')
            axs[x, y].fill_between(df_b['Time'], df_b[column] - df_b[column + '_std'],
                                                       df_b[column] + df_b[column + '_std'],
                                                       color='gray', linewidth=0, alpha=0.3)
            for k, (ip, df_ip) in enumerate(df_delay.groupby('Transmission_To_Human')):
                axs[x, y].plot(df_ip['Time'], df_ip[column], color=cmap(k))
                               # color='xkcd:blue', alpha=(0.1 + 0.9*k/6))

            axs[x, y].xaxis.set_major_locator(majorLocator)
            axs[x, y].xaxis.set_major_formatter(majorFormatter)
            axs[x, y].set_xlim([0, 6 * 365])
            axs[x, y].set_ylim([0.0, 1.0])

            axs[x, y].set_yticks([0.0, 0.5, 1.0])
            axs[x, y].set_yticklabels(['0%', '50%', '100%'])

            axs[x, y].xaxis.set_minor_locator(minorLocator)
            labels = [x for x in range(-1, 7)]

            axs[x, y].xaxis.set_ticklabels(labels)

            if x == 0 and y < 2:
                axs[x, y].xaxis.set_ticklabels([])

            if y > 0:
                axs[x, y].yaxis.set_ticklabels([])

            axs[x, y].set_title('%i%%' % (100-delay*100))

        fig.text(0.01, 0.5, 'True Prevalence', va='center', rotation='vertical', size=24)
        fig.text(0.5, 0.01, 'Year', ha='center', size=24)
        # fig.text(0.85, 0.22, 'Reduction of\ninfectious sporozoites', ha='center', size=16)
        fig.delaxes(axs[1][2])

        cbar_ax = fig.add_axes([.67, .12, .01, .3])
        cbar_ax.yaxis.set_ticks_position('right')

        cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap.reversed(),
                                        norm=norm)
        cb1.set_ticks([4.5 - i for i in range(len(tran_human))])
        cb1.set_ticklabels(['%i%%' % (100 - t*100) for t in tran_human])
        cb1.set_label('Reduction of\ninfectious sporozoites', rotation='horizontal', labelpad=100, va='center')

        plt.subplots_adjust(bottom=0.1, left=0.1)
        plt.suptitle('Average increase in time until sporozoite release')
        plt.savefig(os.path.join(fig_dir, 'prevalence_EIR_%i.pdf' % fig_labels[g]))
        plt.savefig(os.path.join(fig_dir, 'prevalence_EIR_%i.png' % fig_labels[g]))


def plot_allele_frequency():
    df = pd.read_csv(allele_file)
    df = df[df['Baseline'] == 0]
    df = df.groupby(['Larval_Capacity', 'Time'])['a0', 'a1', 'a2', 'a3', 'b0', 'b1', 'b2', 'b3',
    'a0_std', 'a1_std', 'a2_std', 'a3_std', 'b0_std', 'b1_std', 'b2_std', 'b3_std'].apply(np.mean).reset_index()

    lhabitat = list(df['Larval_Capacity'].unique())
    # tran_human = list(df['Transmission_To_Human'].unique())
    EIR = [15, 30, 60]
    fig_labels = dict(zip(lhabitat, EIR))

    columns = ['Driver locus', 'Effector locus']
    genes = {0: 'a', 1: 'b'}

    majorLocator = MultipleLocator(365)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(365 / 12.0)

    colors = ['xkcd:purple', 'xkcd:green', 'xkcd:orange', 'xkcd:gray']
    fig, axs = plt.subplots(3, 2, figsize=(10, 8))

    for i, (g, df_g) in enumerate(df.groupby('Larval_Capacity')):

        for c, column in enumerate(columns):
            x = i
            y = c

            for allele in range(4):
                axs[x, y].plot(df_g['Time'], df_g['%s%i' % (genes[c], allele)], color=colors[allele])
                axs[x, y].fill_between(df_g['Time'],
                                       df_g['%s%i' % (genes[c], allele)] - df_g['%s%i_std' % (genes[c], allele)],
                                       df_g['%s%i' % (genes[c], allele)] + df_g['%s%i_std' % (genes[c], allele)],
                                       color=colors[allele], alpha=0.3)
                    # color='xkcd:blue', alpha=(0.1 + 0.9*k/6))

            axs[x, y].xaxis.set_major_locator(majorLocator)
            axs[x, y].xaxis.set_major_formatter(majorFormatter)
            axs[x, y].set_xlim([0, 6 * 365])
            axs[x, y].set_ylim([0.0, 1.0])

            axs[x, y].set_yticks([0.0, 0.5, 1.0])
            axs[x, y].set_yticklabels(['0%', '50%', '100%'])

            axs[x, y].xaxis.set_minor_locator(minorLocator)
            labels = [x for x in range(-1, 7)]

            axs[x, y].xaxis.set_ticklabels(labels)

            if x < 2:
                axs[x, y].xaxis.set_ticklabels([])

            if y > 0:
                axs[x, y].yaxis.set_ticklabels([])

            if x==0:
                axs[x, y].set_title(columns[c], size=16)

            if y == 1:
                ax2 = axs[x, y].twinx()
                ax2.set_ylabel('Annual EIR\n ~%i' % fig_labels[g], color='b', rotation=0, labelpad=50,
                               fontdict={'size': 16})
                ax2.set_yticks([])

    fig.text(0.01, 0.5, 'Allele Frequencies', va='center', rotation='vertical', size=18)
    fig.text(0.47, 0.11, 'Year', ha='center', size=18)

    driver_alleles = ['Wild type', 'Nuclease or Effector', 'Resistance', 'Loss of gene function']
    custom_lines = [Line2D([0], [0], color=colors[j], lw=2) for j in range(4)]
    fig.legend(custom_lines, [driver_alleles[j] for j in range(4)],
               bbox_to_anchor=(0.75, 0.1), ncol=2, prop={'size': 12})


    plt.subplots_adjust(bottom=0.18, left=0.11, right=0.82)
    plt.savefig(os.path.join(fig_dir, 'allele_freq.pdf'))
    plt.savefig(os.path.join(fig_dir, 'allele_freq.png'))


if __name__ == '__main__':

    plot_prevalence()
    plot_allele_frequency()