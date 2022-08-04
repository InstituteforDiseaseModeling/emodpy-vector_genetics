import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams.update({'font.size': 16})

exp_name = 'sporozoite_reduction'
data_dir = os.path.join(os.path.expanduser('~'), 'Github', 'emodpy-vector_genetics', 'analyzers', 'output')
fig_dir = os.path.join(os.path.expanduser('~'), 'Github', 'emodpy-vector_genetics', 'figures')
os.makedirs(fig_dir, exist_ok=True)

if __name__ == '__main__':

    prev_file = os.path.join(data_dir, 'sporozoite_reduction_prevalence_final.csv')
    allele_file = os.path.join(data_dir, 'sporozoite_reduction_allele_frequency_final.csv')
    incidence_file = os.path.join(data_dir, 'sporozoite_reduction_incidence_final.csv')

    ######## Plot prevalence ########
    df = pd.read_csv(prev_file)
    df_incidence = pd.read_csv(incidence_file)
    df_elimination = df_incidence[df_incidence['Baseline'] == 0]
    df_incidence = df_incidence[df_incidence['Baseline'] == 1]
    df_baseline = df[df['Baseline'] == 1]
    df = df[df['Baseline'] == 0]

    larval_capacity = [7, 7.25, 7.5, 7.75, 8]

    fig, axs = plt.subplots(5, 2, figsize=(8, 10))
    colors = ['xkcd:blue', 'xkcd:orange']
    columns = ['True Prevalence', 'RDT Prevalence']

    for k, l in enumerate(larval_capacity):

        majorLocator = MultipleLocator(365)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(365 / 12.0)

        dftemp = df[df['Larval_Capacity'] == l]
        df_inc = df_incidence[df_incidence['Larval_Capacity'] == l]
        df_elim = df_elimination[df_elimination['Larval_Capacity'] == l]
        incidence_reduction = (1 - df_elim['New Clinical Cases'].values[0]/df_inc['New Clinical Cases'].values[0]) * 100
        df_b = df_baseline[df_baseline['Larval_Capacity'] == l]

        for i, column in enumerate(columns):
            axs[k, i].plot(df_b['Time'], df_b[column], label=column, color='gray')
            axs[k, i].fill_between(df_b['Time'],
                                df_b[column] - df_b[column + '_std'],
                                df_b[column] + df_b[column + '_std'],
                                color='gray', linewidth=0,
                                alpha=0.3)
            axs[k, i].plot(dftemp['Time'], dftemp[column], label=column, color=colors[i])
            axs[k, i].fill_between(dftemp['Time'],
                            dftemp[column] - dftemp[column + '_std'],
                            dftemp[column] + dftemp[column + '_std'],
                            color=colors[i], linewidth=0,
                            alpha=0.3)
            axs[k, i].xaxis.set_major_locator(majorLocator)
            axs[k, i].xaxis.set_major_formatter(majorFormatter)
            axs[k, i].set_xlim([0, 6 * 365])
            axs[k, i].set_ylim([0.0, 1.0])

            if i == 1:
                axs[k, i].text(1.03, 0.3, 'Annual EIR = %i,\nProbability of\nelimination = %i%%,\nIncidence\nreuction (%%) = %i%%'
                               % (np.round(df_inc['Annual EIR'].values[0]), df_elim['Elimination'] * 100,
                                  incidence_reduction), transform=axs[k, i].transAxes, size=12)

            axs[k, i].set_yticks([0.0, 0.5, 1.0])
            axs[k, i].set_yticklabels(['0%', '50%', '100%'])

            axs[k, i].xaxis.set_minor_locator(minorLocator)
            labels = [x for x in range(-1, 7)]

            if k == 4:
                axs[k, i].xaxis.set_ticklabels(labels)
            else:
                axs[k, i].xaxis.set_ticklabels([])

            if i == 1:
                axs[k, i].yaxis.set_ticklabels([])
            if k == 0:
                axs[k, i].set_title(column)

    plt.subplots_adjust(bottom=0.15, right=0.75)

    custom_lines_2 = [Patch(color='gray', alpha=0.3)]
    fig.legend(custom_lines_2, ['Baseline'],
               bbox_to_anchor=(0.55, 0.1), ncol=1)
    # plt.show()
    plt.savefig(os.path.join(fig_dir, 'Prevalence_LC.png'))


    ######## Plot vector genetics ########
    df = pd.read_csv(allele_file)
    df = df[df['Baseline'] == 0]
    df['VectorPopulation_std'] = df['VectorPopulation_std'] / df.groupby(['Baseline', 'Time', 'Larval_Capacity'])[
        'VectorPopulation'].transform('sum')
    df['VectorPopulation'] = df['VectorPopulation'] / df.groupby(['Baseline', 'Time', 'Larval_Capacity'])[
        'VectorPopulation'].transform('sum')

    fig, axs = plt.subplots(5, 2, figsize=(8, 10))
    colors = ['r', 'g', 'b', 'y']
    columns = ['Driver locus', 'Effector locus']
    alleles = {0: 'a', 1: 'b'}

    for k, l in enumerate(larval_capacity):
        df_a = df[df['Larval_Capacity'] == l]
        df_inc = df_incidence[df_incidence['Larval_Capacity'] == l]
        df_elim = df_elimination[df_elimination['Larval_Capacity'] == l]
        incidence_reduction = (1 - df_elim['New Clinical Cases'].values[0] / df_inc['New Clinical Cases'].values[
            0]) * 100

        majorLocator = MultipleLocator(365)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(365 / 12.0)

        for i, column in enumerate(columns):
            for j in range(4):
                dftemp = df_a[df_a['Alleles'] == alleles[i]+str(j)]
                axs[k, i].plot(dftemp['Time'], dftemp['VectorPopulation']*2, label=column, color=colors[j])
                axs[k, i].fill_between(dftemp['Time'],
                                    (dftemp['VectorPopulation'] - dftemp['VectorPopulation_std']) * 2,
                                    (dftemp['VectorPopulation'] + dftemp['VectorPopulation_std']) * 2,
                                    color=colors[j], linewidth=0,
                                    alpha=0.3)
            axs[k, i].xaxis.set_major_locator(majorLocator)
            axs[k, i].xaxis.set_major_formatter(majorFormatter)
            axs[k, i].set_xlim([0, 6 * 365])
            axs[k, i].set_ylim([0.0, 1.0])

            if i == 1:
                axs[k, i].text(1.03, 0.3, 'Annual EIR = %i,\nProbability of\nelimination = %i%%,\nIncidence\nreuction (%%) = %i%%'
                               % (np.round(df_inc['Annual EIR'].values[0]), df_elim['Elimination'] * 100,
                                  incidence_reduction), transform=axs[k, i].transAxes, size=12)

            axs[k, i].set_yticks([0.0, 0.5, 1.0])
            axs[k, i].set_yticklabels(['0%', '50%', '100%'])

            axs[k, i].xaxis.set_minor_locator(minorLocator)
            labels = [x for x in range(-1, 7)]

            if k == 4:
                axs[k, i].xaxis.set_ticklabels(labels)
            else:
                axs[k, i].xaxis.set_ticklabels([])

            if i == 1:
                axs[k, i].yaxis.set_ticklabels([])
            if k == 0:
                axs[k, i].set_title(column)

        # plt.suptitle('Annual EIR = %i, Probability of elimination = %i%%' % (np.round(df_inc['Annual EIR'].values[0]),
        #                                                                  df_elim['Elimination'] * 100))
    plt.subplots_adjust(bottom=0.15, right=0.75)

    driver_alleles = ['Wild type', 'Nuclease or Effector', 'Resistance', 'Loss of gene function']
    custom_lines = [Line2D([0], [0], color=colors[j], lw=2) for j in range(4)]
    fig.legend(custom_lines, [driver_alleles[j] for j in range(4)],
               bbox_to_anchor=(0.7, 0.1), ncol=2, prop={'size': 10})

    # effector_alleles = ['Wild type', 'Effector', 'Resistance', 'Loss of gene function']
    # custom_lines = [Line2D([0], [0], color=colors[j], lw=2) for j in range(4)]
    # fig.legend(custom_lines, [effector_alleles[j] for j in range(4)],
    #            bbox_to_anchor=(0.76, 0.1), ncol=2, prop={'size': 10})
    # plt.show()
    plt.savefig(os.path.join(fig_dir, 'Allele_frequency_LC.png'))
