import os
from functools import \
    partial  # for setting Run_Number. In Jonathan Future World, Run_Number is set by dtk_pre_proc based on generic param_sweep_value...
import sporozoite_delay.manifest as manifest

import emod_api.demographics.Demographics as Demographics

from emodpy_malaria import config as malconf
from emodpy_malaria import vector_config as vecconf
from emodpy_malaria.interventions.treatment_seeking import add_healthseeking
from emodpy_malaria.interventions.mosquitorelease import MosquitoRelease

import emod_api.campaign as camp


def update_sim_bic(simulation, value):
    simulation.task.config.parameters.Base_Infectivity_Constant = value * 0.1
    return {"Base_Infectivity": value}


def update_sim_random_seed(simulation, value):
    simulation.task.config.parameters.Run_Number = value
    return {"Run_Number": value}


def update_serialize(simulation, run_number, serialization=0, sim_duration=40 * 365, serialized_population_path=None):
    if serialization:
        simulation.task.config.parameters.Simulation_Duration = sim_duration
        simulation.task.config.parameters.Serialization_Time_Steps = [sim_duration]
        simulation.task.config.parameters.Serialized_Population_Reading_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Writing_Type = 'TIMESTEP'
        update_sim_random_seed(simulation, value=0)
    else:
        simulation.task.config.parameters.Simulation_Duration = sim_duration
        simulation.task.config.parameters.Serialization_Mask_Node_Read = 0
        simulation.task.config.parameters.Serialization_Mask_Node_Write = 0
        simulation.task.config.parameters.Serialization_Precision = 'REDUCED'
        # simulation.task.config.parameters.Serialization_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Path = os.path.join(serialized_population_path,
                                                                                    'output')
        simulation.task.config.parameters.Serialized_Population_Reading_Type = 'READ'
        simulation.task.config.parameters.Serialized_Population_Writing_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Filenames = ['state-%i.dtk' % sim_duration]
        update_sim_random_seed(simulation, value=run_number)

    return {"Serialization": serialization}


def update_camp_type(simulation, baseline, serialize=0, sim_duration=40 * 365):
    # simulation.task.config.parameters.Run_Number = value
    build_camp_partial = partial(build_camp, serialize=serialize, sim_duration=sim_duration, baseline=baseline)
    simulation.task.create_campaign_from_callback(build_camp_partial)

    return {"Baseline": baseline}


def set_param_fn(config):
    """
    This function is a callback that is passed to emod-api.config to set parameters The Right Way.
    """
    config = malconf.set_team_defaults(config, manifest)
    # config = set_config.set_config(config)

    config.parameters.Base_Rainfall = 150
    config.parameters.Climate_Model = "CLIMATE_CONSTANT"
    config.parameters.Enable_Disease_Mortality = 0
    config.parameters.Enable_Vector_Species_Report = 0
    config.parameters.Demographics_Filenames = ['single_node_demographics.json']
    config.parameters.pop("Serialized_Population_Filenames")

    # # Create new MalariaDrugParams
    # config = set_mdp(config, manifest)
    insecticides = [
        {
            "Name": "pyrethroid",
            "Resistances": []
        }]
    config["parameters"]["Insecticides"] = insecticides

    # Gene drive - using full dictionary here but hope to change to something more elegant in the future
    drivers = {"gambiae": [
        {
            "Alleles_Driven": [
                {
                    "Allele_To_Copy": "a1",
                    "Allele_To_Replace": "a0",
                    "Copy_To_Likelihood": [
                        {
                            "Copy_To_Allele": "a1",
                            "Likelihood": 0.9
                        },
                        {
                            "Copy_To_Allele": "a0",
                            "Likelihood": 0.1
                        }
                    ]
                },
                {
                    "Allele_To_Copy": "b1",
                    "Allele_To_Replace": "b0",
                    "Copy_To_Likelihood": [
                        {
                            "Copy_To_Allele": "b1",
                            "Likelihood": 0.9
                        },
                        {
                            "Copy_To_Allele": "b0",
                            "Likelihood": 0.1
                        }
                    ]
                }
            ],
            "Driver_Type": "INTEGRAL_AUTONOMOUS",
            "Driving_Allele": "a1"
        },
    ]}

    # Vector Genetics
    malconf.set_species(config, ['gambiae'])

    vecconf.set_species_param(config, 'gambiae', 'Larval_Habitat_Types', [{ "Vector_Habitat_Type" : "LINEAR_SPLINE",
        "Capacity_Distribution_Over_Time": {
            "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083,
                      182.5, 212.917, 243.333, 273.75, 304.167, 334.583],
            "Values": [3, 0.8, 1.25, 0.1, 2.7, 10, 6, 35, 2.8, 1.5, 1.6, 2.1]
        },
        "Capacity_Distribution_Number_Of_Years": 1,
        "Max_Larval_Capacity": pow(10, 8)
    }])

    vecconf.add_alleles(['a0', 'a1', 'a2'], [1.0, 0.0, 0.0])
    vecconf.add_alleles(['b0', 'b1', 'b2'], [1.0, 0.0, 0.0])

    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b1"]], "INFECTED_PROGRESS", 0.5)
    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b0"]], "INFECTED_PROGRESS", 0.75)
    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b2"]], "INFECTED_PROGRESS", 0.75)
    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b1"]], "TRANSMISSION_TO_HUMAN", 0.6)
    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b0"]], "TRANSMISSION_TO_HUMAN", 0.8)
    vecconf.add_trait(manifest, [["X", "X"], ["b1", "b2"]], "TRANSMISSION_TO_HUMAN", 0.8)

    vecconf.set_species_drivers(config, drivers)
    vecconf.set_genetics(vecconf.get_species_params(config, "gambiae"), manifest)

    # # Crete new Vector Species Params
    # config = set_vsp(config, manifest)
    return config


def build_camp(serialize=0, sim_duration=40 * 365, baseline=1):
    """
    Build a campaign input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    """

    # This isn't desirable. Need to think about right way to provide schema (once)
    camp.schema_path = manifest.schema_file

    # print( f"Telling emod-api to use {manifest.schema_file} as schema." )
    # camp.add(bednet.Bednet(camp, start_day=start_day_in, coverage=0.5, killing_eff=0.5, blocking_eff=0.5, usage_eff=0.5,
    #                        insecticide="pyrethroid"), first=True)
    if not serialize:
        add_healthseeking(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 0.8, "agemin": 15,
                                          "agemax": 70, "seek": 0.4, "rate": 0.3},
                                         {"trigger": "NewSevereCase", "coverage": 0.8, "seek": 0.6, "rate": 0.5}],
                          drug=['Artemether', 'Lumefantrine'],
                          start_day=0,
                          broadcast_event_name='Received_Treatment'
                          )
        if not baseline:
            camp.add(MosquitoRelease(camp, start_day=180, species='gambiae', genome=[["X", "X"], ["a1", "a1"],
                                                                                     ["b1", "b1"]], number=100))

    else:
        add_healthseeking(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 0.8, "agemin": 15,
                                          "agemax": 70, "seek": 0.4, "rate": 0.3},
                                         {"trigger": "NewSevereCase", "coverage": 0.8, "seek": 0.6, "rate": 0.5}],
                          drug=['Artemether', 'Lumefantrine'],
                          start_day=sim_duration - 10 * 365,
                          broadcast_event_name='Received_Treatment'
                          )
    return camp


def build_demog():
    """
    Build a demographics input file for the DTK using emod_api.
    Right now this function creates the file and returns the filename. If calling code just needs an asset that's fine.
    Also right now this function takes care of the config updates that are required as a result of specific demog settings. We do NOT want the emodpy-disease developers to have to know that. It needs to be done automatically in emod-api as much as possible.
    TBD: Pass the config (or a 'pointer' thereto) to the demog functions or to the demog class/module.

    """
    demog = Demographics.from_file("C:\\Users\\pselvaraj\\Github\\emodpy-vector_genetics\\input_files\\single_node_demographics.json")

    return demog


def rvg_config_builder(params):
    params.Combine_Similar_Genomes = True
    params.Include_Vector_State_Columns = False
    params.Species = 'gambiae'
    # params.Alleles_For_Stratification = [ "tom", "dick" ] # this works
    # params.Specific_Genome_Combinations_For_Stratification =
    """
    E.g.,
    [
        {
            "Allele_Combination": [
                ["X", "X"],
                ["a1", "*"]
            ]
        },
        {
            "Allele_Combination": [
                ["X", "X"],
                ["a0", "a0"]
            ]
        }
    ]
    """
    return params
