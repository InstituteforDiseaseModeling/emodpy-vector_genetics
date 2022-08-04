import os
import numpy as np
import pandas as pd
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser

column_types = dict(Time='uint16',
                    time='uint16',
                    Run_Number='uint8',
                    Node='uint8',
                    node='uint8',
                    NodeID='uint8'
                    )


def optimize_dataframe(simdata):
    """
    Set the data to the optimal size for their types. For examples, we know this experiment does not have
    more than 255 nodes so we should assign the NodeId columns to an uint8 since it takes less memory than the
    default uint64
    :param simdata:
    :return:
    """
    for column, item_type in column_types.items():
        if column in simdata.columns:
            simdata[column] = simdata[column].astype(item_type)
    return simdata


class AlleleFreqAnalyzer(BaseAnalyzer):

    def __init__(self, exp_name, sweep_variables, working_dir='.'):
        super(AlleleFreqAnalyzer, self).__init__(
            working_dir=working_dir,
            # filenames=['output/ReportVectorGenetics_gambiae_ALLELE_FREQ.csv']
            filenames=['output/ReportVectorGenetics_gambiae_Female_ALLELE_FREQ.csv']
        )
        self.exp_name = exp_name
        self.sweep_variables = sweep_variables
        self.restrict_nodes = 1
        if self.restrict_nodes == 1:
            # self.desired_nodes = [3, 77, 112, 144]
            self.desired_nodes = [1]
        self.num_loci = 3  # including sex chromosome
        if self.num_loci == 2:
            self.channels = ['a0', 'a1', 'a2',
                             'Total Female Population']
        elif self.num_loci == 3:
            self.channels = ['a0', 'a1', 'a2', 'a3', 'b0', 'b1', 'b2', 'b3',
                             'Total Female Population']
        elif self.num_loci == 4:
            self.channels = ['a0', 'a1', 'a2', 'a3',
                             'b0', 'b1', 'b2', 'b3',
                             'c0', 'c1', 'c2', 'c3',
                             'Total Female Population']
        elif self.num_loci == 5:
            self.channels = ['a0', 'a1', 'a2', 'a3',
                             'b0', 'b1', 'b2', 'b3',
                             'c0', 'c1', 'c2', 'c3',
                             'd0', 'd1', 'd2', 'd3',
                             'Total Female Population']
        self.reporter_sex = 'female'  # choose: female, male, both
        self.output_type = 'meanstd'  # choose: sims, meanstd
        if self.output_type == 'sims':
            self.output_fname = os.path.join(self.working_dir, "%s_spatial_avg_allele_freqs_indiv_sims" % self.exp_name)
        elif self.output_type == 'meanstd':
            self.output_fname = os.path.join(self.working_dir, "%s_avg_allele_freqs" % self.exp_name)

    def select_simulation_data(self, data, simulation):
        simdata = []

        if self.reporter_sex == 'female':
            datatemp = data['output/ReportVectorGenetics_gambiae_Female_ALLELE_FREQ.csv']
            datatemp = datatemp[datatemp['Alleles'] != 'Y']
            if self.restrict_nodes == 1:
                datatemp = datatemp[datatemp['NodeID'].isin(self.desired_nodes)]
            alleles = datatemp['Alleles'].unique()
            alleles = alleles[alleles != 'X']
            datatemp = datatemp.pivot_table('VectorPopulation', ['Time', 'NodeID'], 'Alleles').reset_index()
            datatemp = datatemp.groupby('Time').sum().reset_index()
            pop_col_name = 'Total Female Population'
            datatemp[pop_col_name] = datatemp['X'] / 2
            datatemp.drop(columns=['NodeID', 'X'], inplace=True)
            for allele in alleles:
                datatemp[allele] = datatemp[allele] / (2 * datatemp[pop_col_name])

        elif self.reporter_sex == 'male':
            datatemp = data['output/ReportVectorGenetics_gambiae_Male_ALLELE_FREQ.csv']
            datatemp = datatemp[datatemp['Alleles'] != 'Y']
            if self.restrict_nodes == 1:
                datatemp = datatemp[datatemp['NodeID'].isin(self.desired_nodes)]
            alleles = datatemp['Alleles'].unique()
            alleles = alleles[alleles != 'X']
            datatemp = datatemp.pivot_table('STATE_MALE', ['Time', 'NodeID'],
                                            'Alleles').reset_index()  # NOT SURE ABOUT STATE_MALE HERE
            datatemp = datatemp.groupby('Time').sum().reset_index()
            pop_col_name = 'Total Male Population'
            datatemp[pop_col_name] = datatemp['X']
            datatemp.drop(columns=['NodeID', 'X'], inplace=True)
            for allele in alleles:
                datatemp[allele] = datatemp[allele] / (2 * datatemp[pop_col_name])

        elif self.reporter_sex == 'both':
            datatemp = data['output/ReportVectorGenetics_gambiae_ALLELE_FREQ.csv']
            if self.restrict_nodes == 1:
                datatemp = datatemp[datatemp['NodeID'].isin(self.desired_nodes)]
            alleles = datatemp['Alleles'].unique()
            alleles = alleles[(alleles != 'X') & (alleles != 'Y')]
            datatemp = datatemp.pivot_table('VectorPopulation', ['Time', 'NodeID'], 'Alleles').reset_index()
            datatemp = datatemp.groupby('Time').sum().reset_index()
            pop_col_name = 'Total Female Population'
            datatemp[pop_col_name] = datatemp['X'] / 2
            datatemp.drop(columns=['NodeID', 'X'], inplace=True)
            for allele in alleles:
                datatemp[allele] = datatemp[allele] / (2 * datatemp[pop_col_name])
            # datatemp = datatemp.pivot_table(['VectorPopulation', 'STATE_MALE'], ['Time', 'NodeID'], 'Alleles').reset_index()
            # datatemp = datatemp.groupby('Time').sum().reset_index()
            # datatemp['Total Male Population'] = datatemp['STATE_MALE']['Y']
            # datatemp['Total Female Population'] = datatemp['VectorPopulation']['X'] / 2
            # for allele in alleles:
            #     datatemp['male ' + allele] = datatemp['STATE_MALE'][allele] / (2 * datatemp['Total Male Population'])
            # for allele in alleles:
            #     datatemp['female ' + allele] = datatemp['VectorPopulation'][allele] / (2 * datatemp['Total Female Population'])
            # datatemp = datatemp.drop('STATE_MALE', axis=1, level=0).drop(
            #     'VectorPopulation', axis=1, level=0).drop(
            #     'NodeID', axis=1, level=0)

        simdata.append(datatemp)
        simdata = pd.concat(simdata)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
            else:
                simdata[sweep_var] = 0

        simdata.reset_index(drop=True)
        return optimize_dataframe(simdata)

    def finalize(self, d: pd.DataFrame, partition_vars: list, partition_vars_vals: list):

        fname_suffix = ''
        for ipv, partition_var in enumerate(partition_vars):
            fname_suffix = fname_suffix + '_' + partition_var + str(partition_vars_vals[ipv])

        if self.output_type == 'sims':
            d.to_csv(self.output_fname + fname_suffix + '.csv', index=False)

        elif self.output_type == 'meanstd':
            d_mean = d.groupby(self.sweep_variables + ['Time'])[self.channels].apply(
                np.mean).reset_index()
            d_std = d.groupby(self.sweep_variables + ['Time'])[self.channels].apply(
                np.std).reset_index()
            for channel in self.channels:
                d_mean[channel + '_std'] = d_std[channel]
            d_mean.to_csv(self.output_fname + fname_suffix + '.csv', index=False)


if __name__ == "__main__":

    SetupParser.default_block = 'HPC'
    # SetupParser.default_block = 'LOCAL'
    SetupParser.init()

    # projectdir = os.path.join(os.path.expanduser('~'), 'Dropbox (IDM)', 'Malaria Team Folder', 'projects',
    #                           'Vector_genetics',
    #                           'Data', 'simulation_data', 'Establishment_Burkina')
    # out_dir = os.path.join(projectdir)
    out_dir = 'C:\\Users\\pselvaraj\\Github\\emodpy - vector_genetics\\analyzers\output'

    # experiments = {
    #     "spatial_vector_genetics_classicdrive3allele_VC_and_GM_sweep_vmiginoutfrac":
    #         "b1db125b-cf53-eb11-a2dd-c4346bcb7271"
    # }
    experiments = {
        "spatialinside_vector_genetics_classic3allele_VC_and_GM_vmignegexpovnnb05_sweep_rc_d_rr0_sne":
            "cbc92044-45be-eb11-a9ec-b88303911bc1"
    }

    # sweep_vars = ['vmigin_frac', 'vmigout_frac']
    sweep_vars = ['Baseline', 'Run_Number', 'Larval_Capacity', 'Transmission_To_Human', 'Infected_Progress']

    for expt_name, exp_id in experiments.items():
        am = AnalyzeManager(exp_list=exp_id,
                            analyzers=[
                                AlleleFreqAnalyzer(exp_name='spatialinside_vector_genetics_classic3allele_VC_and_GM_vmignegexpovnnb05_sweep_rc_d_rr0_sne',
                                                   working_dir=out_dir,
                                                   sweep_variables=sweep_vars
                                                   )
                            ],
                            force_analyze=False)

        print(am.experiments)
        am.analyze()