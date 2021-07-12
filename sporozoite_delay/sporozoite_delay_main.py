#!/usr/bin/env python3

import pathlib # for a join

# idmtools ...
from idmtools.builders import SimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from idmtools.utils.filter_simulations import FilterItem

# emodpy
from emodpy.emod_task import EMODTask
from emodpy_malaria.reporters.builtin import ReportVectorGenetics

from sporozoite_delay.helpers import *
import sporozoite_delay.params as params
import sporozoite_delay.manifest as manifest

from functools import partial
import pandas as pd
from uuid import UUID


def get_serialization_paths(platform, serialization_exp_id):
    exp = Experiment.from_id(serialization_exp_id, children=False)
    exp.simulations = platform.get_children(exp.id, exp.item_type,
                                            children=["tags", "configuration", "files", "hpc_jobs"])

    sim_dict = {'Larval_Capacity': [], 'Outpath': []}
    for simulation in exp.simulations:
        # if simulation.tags['Run_Number'] == 0:
        string = simulation.get_platform_object().hpc_jobs[0].working_directory.replace('internal.idm.ctr', 'mnt')
        string = string.replace('\\', '/')
        string = string.replace('IDM2', 'idm2')

        sim_dict['Larval_Capacity'] += [float(simulation.tags['Larval_Capacity'])]
        sim_dict['Outpath'] += [string]

    df = pd.DataFrame(sim_dict)
    return df


def general_sim(serialization=0, serialized_exp_id=None):
    """
    This function is designed to be a parameterized version of the sequence of things we do 
    every time we run an emod experiment. 
    """

    # Create a platform
    # Show how to dynamically set priority and node_group
    platform = Platform("SLURM")

    # create EMODTask 
    print("Creating EMODTask (from files)...")
    
    task = EMODTask.from_default2(
            config_path="my_config.json",
            eradication_path=manifest.eradication_path,
            ep4_custom_cb=None,
            campaign_builder=None,
            schema_path=manifest.schema_file,
            param_custom_cb=set_param_fn,
            demog_builder=None,
        )

    # Create simulation sweep with builder
    builder = SimulationBuilder()

    # Add asset
    task.common_assets.add_asset("C:\\Users\\pselvaraj\\Github\\emodpy-vector_genetics\\input_files\\single_node_demographics.json")

    if serialized_exp_id:
        serialized_population_path_df = get_serialization_paths(platform=platform, serialization_exp_id=serialized_exp_id)

        func = partial(update_serialize, serialization=serialization, sim_duration=10 * 365,
                       serialized_population_path_df=serialized_population_path_df)
        builder.add_sweep_definition(func, [7.0, 7.25, 7.5, 7.75, 8.0])

        builder.add_sweep_definition(update_sim_random_seed, range(params.nSims))
        func = partial(update_camp_type, serialize=serialization, sim_duration=10 * 365)
        builder.add_sweep_definition(func, [True, False])
        exp_name = params.exp_name

        # Add reporter
        reporter = ReportVectorGenetics()  # Create the reporter
        reporter.config(rvg_config_builder, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)  # Add thre reporter

    else:
        func = partial(update_serialize, serialization=serialization, sim_duration=40 * 365,
                       serialized_population_path_df=None)
        builder.add_sweep_definition(func, [7.0, 7.25, 7.5, 7.75, 8.0])
        builder.add_sweep_definition(update_sim_random_seed, [0])
        func = partial(update_camp_type, serialize=serialization, sim_duration=40 * 365)
        builder.add_sweep_definition(func, [1])
        exp_name = params.exp_name + '_serialization'

    # create experiment from builder
    print( f"Prompting for COMPS creds if necessary..." )
    experiment = Experiment.from_builder(builder, task, name=exp_name)

    # The last step is to call run() on the ExperimentManager to run the simulations.
    experiment.run(wait_until_done=True, platform=platform)

    # Check result
    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.\n")
        exit()

    print(f"Experiment {experiment.uid} succeeded.")

    # Save experiment id to file
    with open("COMPS_ID", "w") as fd:
        fd.write(experiment.uid.hex)
    print()
    print(experiment.uid.hex)


if __name__ == "__main__":
    # TBD: user should be allowed to specify (override default) erad_path and input_path from command line 
    # plan = EradicationBambooBuilds.MALARIA_LINUX
    # print("Retrieving Eradication and schema.json from Bamboo...")
    # get_model_files( plan, manifest )
    # print("...done.")

    serialization = 0
    serialization_experiment_id = '9343ae74-e1de-eb11-a9ec-b88303911bc1'
    # serialization = 1
    # serialization_experiment_id = None
    general_sim(serialization=serialization, serialized_exp_id=serialization_experiment_id)
