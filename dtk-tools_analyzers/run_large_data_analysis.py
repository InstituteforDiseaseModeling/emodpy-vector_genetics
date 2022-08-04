import os
import sys

sys.path.append(os.path.dirname(__file__))

from partitioning_analyzer import PartitioningDataAnalyzeManager

from allele_frequency_analyzer import AlleleFreqAnalyzer


if __name__ == "__main__":

    experiments = {
        "sporozoite_delay_gene_drive_LC_sweep":
            "8bbe4c48-6326-ec11-9ecd-9440c9bee941"
    }

    sweep_vars = ['Baseline', 'Larval_Capacity', 'Transmission_To_Human', 'Infected_Progress']

    for expt_name, exp_id in experiments.items():
        am = PartitioningDataAnalyzeManager(exp_list=exp_id,
                                            partitionable_columns=sweep_vars,
                                            analyzers=
                                            [
                                                AlleleFreqAnalyzer(
                                                    exp_name=expt_name,
                                                    sweep_variables=sweep_vars
                                                ),
                                            ]
                                            )
        print(am.experiments)
        am.analyze()